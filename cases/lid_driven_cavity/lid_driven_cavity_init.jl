## Loading packages.
using MicroHH
using Tullio
using Statistics


## Loading settings.
include("lid_driven_cavity_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("lid_driven_cavity", n_domains, settings, float_type)


## Create the initials fields.
f = m.fields[1]; g = m.grid[1]
u = @view f.u[g.is+1:g.ie, g.js:g.je, g.ks:g.ke]
w = @view f.w[g.is:g.ie, g.js:g.je, g.ks+1:g.ke]
u .+= 1e-4 .* rand(float_type, size(u))
w .+= 1e-4 .* rand(float_type, size(w))
u .-= mean(u, dims=[1, 2])
w .-= mean(w, dims=[1, 2])

u_bot = @view f.u_bot[g.is:g.ie, g.js:g.je]
u_top = @view f.u_top[g.is:g.ie, g.js:g.je]
v_bot = @view f.v_bot[g.is:g.ie, g.js:g.je]
v_top = @view f.v_top[g.is:g.ie, g.js:g.je]
f.u_bot[:, :] .= 0.
# f.u_top[:, :] .= 1.
f.u_top[:, :] .= 0.
f.v_bot[:, :] .= 0.
f.v_top[:, :] .= 0.


## Save the restart files.
save_model(m)

