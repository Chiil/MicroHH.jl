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
fig = Figure(resolution=(900, 800))
f1 = m.fields[1]; g1 = m.grid[1]
x1 = @view g1.x[g1.is:g1.ie]
y1 = @view g1.y[g1.js:g1.je]
z1 = @view g1.z[g1.js:g1.je]
s_bot1 = @view f1.s_bot[g1.is:g1.ie, g1.js:g1.je]
s1 = @view f1.s[g1.is:g1.ie, g1.js:g1.je, g1.ks:g1.ke]
node11 = Node(s_bot1 .- mean(s_bot1))
node12 = Node(s1[:, 1, :] .- mean(s1, dims=(1, 2))[:, 1, :])
h11 = heatmap(fig[1, 1], x1, y1, node11, colorrange=(-0.4, 0.4))
h12 = heatmap(fig[1, 2], x1, z1, node12, colorrange=(-0.4, 0.4))

f2 = m.fields[2]; g2 = m.grid[2]
x2 = @view g2.x[g2.is:g2.ie]
y2 = @view g2.y[g2.js:g2.je]
z2 = @view g2.z[g2.js:g2.je]
s_bot2 = @view f2.s_bot[g2.is:g2.ie, g2.js:g2.je]
s2 = @view f2.s[g2.is:g2.ie, g2.js:g2.je, g2.ks:g2.ke]
node21 = Node(s_bot2 .- mean(s_bot2))
node22 = Node(s2[:, 1, :] .- mean(s2, dims=(1, 2))[:, 1, :])
h21 = heatmap(fig[2, 1], x2, y2, node21, colorrange=(-0.4, 0.4))
h22 = heatmap(fig[2, 2], x2, z2, node22, colorrange=(-0.4, 0.4))
fig


## Run the model.
while in_progress
    global in_progress = step_model!(m)
    node11[] = s_bot1 .- mean(s_bot1)
    node21[] = s_bot2 .- mean(s_bot2)
    node12[] = s1[:, 1, :] .- mean(s1, dims=(1, 2))[:, 1, :]
    node22[] = s2[:, 1, :] .- mean(s2, dims=(1, 2))[:, 1, :]
end
