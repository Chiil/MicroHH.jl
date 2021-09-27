## Loading packages.
using MicroHH
using Tullio
using Statistics
using PyPlot
pygui(:qt5)


## Solve the Taylor Green vortex
function solve_taylorgreen(itot)
    # Loading settings.
    d_grid = Dict{String, Any}()
    d_grid["itot"] = itot
    d_grid["jtot"] = 1
    d_grid["ktot"] = itot ÷ 2

    d_grid["xsize"] = 1.
    d_grid["ysize"] = 1.
    d_grid["zsize"] = 0.5

    zsize = d_grid["zsize"]; ktot = d_grid["ktot"]
    dz = zsize / ktot
    d_grid["z"] = range(0.5*dz, step=dz, length=ktot) |> collect

    d_fields = Dict{String, Any}()
    d_fields["visc"] = (8 * pi^2 * 100)^(-1)
    d_fields["alpha"] = 0.

    d_boundary = Dict{String, Any}()
    d_boundary["mom_bot_type"] = "Neumann"
    d_boundary["mom_top_type"] = "Neumann"
    d_boundary["s_bot_type"] = "Neumann"
    d_boundary["s_top_type"] = "Neumann"

    d_timeloop = Dict{String, Any}()
    d_timeloop["start_time"] = 0.
    d_timeloop["end_time"] = 1.
    d_timeloop["save_time"] = 100.
    d_timeloop["check_time"] = 0.1
    d_timeloop["dt"] = 0.005 * (128 / itot)

    d_multidomain = Dict{String, Any}()
    d_multidomain["enable_nudge"] = false

    settings = [ Dict(
        "grid" => d_grid,
        "fields" => d_fields,
        "boundary" => d_boundary,
        "timeloop" => d_timeloop,
        "multidomain" => d_multidomain) ]


    # Initialize the model.
    n_domains = 1
    m = Model("taylorgreen", n_domains, settings, Float64)


    # Create the initials fields.
    f = m.fields[1]; g = m.grid[1]
    u = @view f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke]
    w = @view f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh]
    p = @view f.p[g.is:g.ie, g.js:g.je, g.ks:g.ke]
    x = @view g.x[g.is:g.ie]
    xh = @view g.xh[g.is:g.ie]
    z = @view g.z[g.ks:g.ke]
    zh = @view g.zh[g.ks:g.keh]
    @tullio u[i, j, k] = sin(2pi*xh[i])*cos(2pi*z[k])
    @tullio w[i, j, k] = -cos(2pi*x[i])*sin(2pi*zh[k])

    f_t = exp(-8. * pi^2 * d_fields["visc"] * d_timeloop["end_time"])
    u_ref = u[:, :, :] * f_t
    w_ref = w[:, :, :] * f_t
    p_ref = zeros(size(u_ref))
    @tullio p_ref[i, j, k] = (1/4 * (cos(4pi*x[i]) + cos(4pi*z[k])) - 1/4) * f_t^2


    # Run the model
    in_progress = prepare_model!(m)
    while in_progress
        in_progress = step_model!(m)
    end

    dx = 1. / itot
    u_err = sum(dx^2 * abs.(u .- u_ref))
    w_err = sum(dx^2 * abs.(w .- w_ref))
    p_err = sum(dx^2 * abs.(p .- p_ref))

    return u_err, w_err, p_err
end


## Solve the case at multiple resolutions.
u016, w016, p016 = solve_taylorgreen(16)
u032, w032, p032 = solve_taylorgreen(32)
u064, w064, p064 = solve_taylorgreen(64)
u128, w128, p128 = solve_taylorgreen(128)
u256, w256, p256 = solve_taylorgreen(256)

dxs = [ 1/16, 1/32, 1/64, 1/128, 1/256 ]
u_errs = [ u016, u032, u064, u128, u256 ]
w_errs = [ w016, w032, w064, w128, w256 ]
p_errs = [ p016, p032, p064, p128, p256 ]


## Plot the results
figure()
loglog(dxs, u_errs, "C0o-", label="u")
loglog(dxs, w_errs, "C1o-", label="w")
loglog(dxs, p_errs, "C2o-", label="p")
loglog(dxs, 0.6*dxs.^2, "k:", label="2nd")
legend(loc=0, frameon=false)
xlabel(L"$\Delta$ x")
ylabel("Error")
display(gcf())
show()
