## Loading packages.
using MicroHH
using Tullio
using Statistics


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings, float_type)


## Create the initials fields.
f = m.fields[1]; g = m.grid[1]
s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]
z = @view g.z[g.ks:g.ke]
rand2d = rand(g.imax, g.jmax)
rand2d .-= mean(rand2d)
s[:, :, 1] .+= rand2d[:, :]
f.s_gradbot[:, :] .= - 0.1 / settings[1]["fields"]["visc"]
f.s_gradtop[:, :] .= 0.003
@tullio s[i, j, k] += 0.003 * z[k]


## Save the restart files.
save_model(m)

