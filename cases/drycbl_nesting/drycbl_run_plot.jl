## Loading packages.
using MicroHH
using Tullio
using Statistics
using GLMakie


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 2
m = Model("drycbl", n_domains, settings, float_type)
load_model!(m)
in_progress = prepare_model!(m)


## Make the figure.
fig = Figure(resolution=(800, 500))

f = m.fields[1]; g = m.grid[1]
x = @view g.x[g.is:g.ie]
y = @view g.y[g.js:g.je]
z = @view g.z[g.ks:g.ke]
s_bot = @view f.s_bot[g.is:g.ie, g.js:g.je]
s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]
s_ref = f.s_ref[g.ks:g.ke]

node1 = Observable(s_bot .- mean(s_bot))
node2 = Observable(s[:, 1, :] .- mean(s, dims=(1, 2))[:, 1, :])
h1 = heatmap(fig[1, 1], x, y, node1, colorrange=(-1.0, 1.0))
h2 = heatmap(fig[1, 2], x, z, node2, colorrange=(-1.0, 1.0))

f2 = m.fields[2]; g2 = m.grid[2]
x2 = @view g2.x[g2.is:g2.ie]
y2 = @view g2.y[g2.js:g2.je]
z2 = @view g2.z[g2.ks:g2.ke]
s_bot2 = @view f2.s_bot[g2.is:g2.ie, g2.js:g2.je]
s2 = @view f2.s[g2.is:g2.ie, g2.js:g2.je, g2.ks:g2.ke]
s_ref2 = f2.s_ref[g2.ks:g2.ke]

node3 = Observable(s_bot2 .- mean(s_bot2))
node4 = Observable(s2[:, 1, :] .- mean(s2, dims=(1, 2))[:, 1, :])
h3 = heatmap!(fig[1, 1], x2, y2, node3, colorrange=(-1.0, 1.0))
h4 = heatmap!(fig[1, 2], x2, z2, node4, colorrange=(-1.0, 1.0))

display(fig)


## Run the model.
while in_progress
    global in_progress = step_model!(m)
    node1[] = s_bot .- mean(s_bot)
    node2[] = s[:, 1, :] .- mean(s, dims=(1, 2))[:, 1, :]
    node3[] = s_bot2 .- mean(s_bot2)
    node4[] = s2[:, 1, :] .- mean(s2, dims=(1, 2))[:, 1, :]
end
