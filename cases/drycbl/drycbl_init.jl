## Loading packages.
using MicroHH
using Tullio
using Statistics


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings)


## Create the initials fields.
f = m.fields[1]; g = m.grid[1]
s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]
x = @view g.x[g.is:g.ie]
y = @view g.y[g.js:g.je]
z = @view g.z[g.ks:g.ke]
rand2d = rand(g.imax, g.jmax)
rand2d .-= mean(rand2d)
s[:, :, 1] .+= rand2d[:, :]

s_gradbot = @view f.s_gradbot[g.is:g.ie, g.js:g.je]
s_gradbot[:, :] .= - 0.1 / settings[1]["fields"]["visc"]

f.s_gradtop[:, :] .= 0.003
@tullio s[i, j, k] += 0.003 * z[k]

# Add hot plume.
smoke_scalar = false

if smoke_scalar
    smoke_gradbot = @view f.scalars_gradbot["smoke"][g.is:g.ie, g.js:g.je]
    x_2d = reshape(x, (size(x, 1), 1))
    y_2d = reshape(y, (1, size(y, 1)))
    
    σ = 100
    @. smoke_gradbot += - 0.1 / settings[1]["fields"]["visc"] * exp( -(x_2d - g.xsize/2)^2 / σ^2 ) * exp( -(y_2d - g.ysize/2)^2 / σ^2 )
    @. s_gradbot += - 0.4 / settings[1]["fields"]["visc"] * exp( -(x_2d - g.xsize/2)^2 / σ^2 ) * exp( -(y_2d - g.ysize/2)^2 / σ^2 )
end
    

## Save the restart files.
save_model(m)

