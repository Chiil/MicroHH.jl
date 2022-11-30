## Loading packages.
using MicroHH
using Tullio
using Statistics


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 2
m = Model("drycbl", n_domains, settings, float_type)


## Create the initials fields.
for i in 1:n_domains
    f = m.fields[i]; g = m.grid[i]
    x = @view g.x[g.is:g.ie]
    z = @view g.z[g.ks:g.ke]
    s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]
    s_gradbot = @view f.s_gradbot[g.is:g.ie, g.js:g.je]

    rand2d = 1e-2 .* rand(g.imax, g.jmax)
    rand2d .-= mean(rand2d)
    s[:, :, 1] .+= rand2d[:, :]
    @tullio s[i, j, k] += 0.003 * z[k]

    # Set the reference profile.
    s_ref = @view f.s_ref[g.ks:g.ke]
    s_ref[:] .= 0.003 * z[:]

    # Set the surface fluxes based on an absolute coordinate
    sgrad = - 0.1 / settings[i]["fields"]["visc"]
    xsize_max = settings[1]["grid"]["xsize"]
    @tullio s_gradbot[i, j] = 2 * (x[i] / xsize_max) * sgrad

    f.s_gradtop[:, :] .= 0.003
end


## Save the restart files.
save_model(m)

