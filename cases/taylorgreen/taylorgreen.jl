## Loading packages.
using MicroHH
using Tullio
using Statistics
using PyPlot
pygui(:qt5)


## Loading settings.
d_grid = Dict{String, Any}()
d_grid["itot"] = 128
d_grid["jtot"] = 1
d_grid["ktot"] = 64

d_grid["xsize"] = 1.
d_grid["ysize"] = 1.
d_grid["zsize"] = 0.5

zsize = d_grid["zsize"]; ktot = d_grid["ktot"]
dz = zsize / ktot
d_grid["z"] = range(0.5*dz, step=dz, length=ktot) |> collect

d_fields = Dict{String, Any}()
d_fields["visc"] = (8 * pi^2 * 1000)^(-1)
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
d_timeloop["check_time"] = 0.01
d_timeloop["dt"] = 0.001

settings = [ Dict(
    "grid"     => d_grid,
    "fields"   => d_fields,
    "boundary" => d_boundary,
    "timeloop" => d_timeloop) ]


## Initialize the model.
n_domains = 1
m = Model("taylorgreen", n_domains, settings, Float32)


## Create the initials fields.
f = m.fields[1]; g = m.grid[1]
u = @view f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke]
w = @view f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh]
x = @view g.x[g.is:g.ie]
xh = @view g.xh[g.is:g.ie]
z = @view g.z[g.ks:g.ke]
zh = @view g.zh[g.ks:g.keh]
@tullio u[i, j, k] = sin(2pi*x[i])*cos(2pi*z[k])
@tullio w[i, j, k] = -cos(2pi*x[i])*sin(2pi*zh[k])

u_ref = u[:, :, :] * exp(-8. * pi^2. * d_fields["visc"] * d_timeloop["end_time"])
w_ref = w[:, :, :] * exp(-8. * pi^2. * d_fields["visc"] * d_timeloop["end_time"])


## Run the model
prepare_model!(m)
in_progress = true
while in_progress
    global in_progress = step_model!(m)
end


## Plot the end result.
u_err = u .- u_ref
w_err = w .- w_ref
println(maximum(u_err))
println(maximum(w_err))

figure()
subplot(131)
pcolormesh(xh, z, u[:, 1, :]')
colorbar()
subplot(132)
pcolormesh(xh, z, u_ref[:, 1, :]')
colorbar()
subplot(133)
pcolormesh(xh, z, u_err[:, 1, :]')
colorbar()
tight_layout()
display(gcf())

figure()
subplot(131)
pcolormesh(x, zh, w[:, 1, :]')
colorbar()
subplot(132)
pcolormesh(x, zh, w_ref[:, 1, :]')
colorbar()
subplot(133)
pcolormesh(x, zh, w_err[:, 1, :]')
colorbar()
tight_layout()
display(gcf())
