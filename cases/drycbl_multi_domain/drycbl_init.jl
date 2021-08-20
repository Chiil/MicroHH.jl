## Loading packages.
using MicroHH
using Tullio
using Statistics


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 2
m = Model("drycbl", n_domains, settings)


## Create the initials fields.
for i in 1:n_domains
    f = m.fields[i]; g = m.grid[i]
    s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]
    z = range(g.dz[1]/2, g.zsize, step=g.dz[1]) |> collect
    rand2d = rand(g.itot, g.jtot)
    rand2d .-= mean(rand2d)
    s[:, :, 1] .+= rand2d[:, :]
    f.s_gradbot[:, :] .= - 0.1 / settings[i]["fields"]["visc"]
    @tullio s[i, j, k] += 0.003 * z[k]
end


## Save the restart files.
save_model(m)

