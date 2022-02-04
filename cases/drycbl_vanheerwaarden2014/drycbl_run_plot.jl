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
z = @view g.z[g.ks:g.ke]
zh = @view g.zh[g.ks:g.keh]
w = @view f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh]
s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]

fig = Figure(resolution=(800, 500))
node1 = Observable(w[:, 1, :])
h1 = heatmap(fig[1, 1], x, zh, node1, colorrange=(-0.2, 0.2))
node2 = Observable(s[:, 1, :] .- mean(s, dims=(1, 2))[:, 1, :])
h2 = heatmap(fig[1, 2], x, z, node2, colorrange=(-0.3, 0.3))
fig


## Run the model.
while in_progress
    global in_progress = step_model!(m)
    node1[] = w[:, 1, :]
    node2[] = s[:, 1, :] .- mean(s, dims=(1, 2))[:, 1, :]
end
