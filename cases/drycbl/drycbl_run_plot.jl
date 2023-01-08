## Loading packages.
using MicroHH
using Tullio
using Statistics
using GLMakie


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings, float_type)
load_model!(m)
in_progress = prepare_model!(m)


## Make the figure.
f = m.fields[1]; g = m.grid[1]
x = @view g.x[g.is:g.ie]
y = @view g.y[g.js:g.je]
z = @view g.z[g.ks:g.ke]
s_bot = @view f.s_bot[g.is:g.ie, g.js:g.je]
s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]

fig = Figure(resolution=(800, 500))
node1 = Observable(s_bot .- mean(s_bot))
h1 = heatmap(fig[1, 1], x, y, node1, colorrange=(-0.4, 0.4))
node2 = Observable(s[:, 1, :] .- mean(s, dims=(1, 2))[:, 1, :])
h2 = heatmap(fig[1, 2], x, z, node2, colorrange=(-0.4, 0.4))
display(fig)


## Run the model.
while in_progress
    global in_progress = step_model!(m)
    node1[] = s_bot .- mean(s_bot)
    node2[] = s[:, 1, :] .- mean(s, dims=(1, 2))[:, 1, :]
end
