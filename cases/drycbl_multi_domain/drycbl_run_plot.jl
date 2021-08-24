## Loading packages.
using MicroHH
using Tullio
using Statistics
using GLMakie


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 2
m = Model("drycbl", n_domains, settings, Float32)
load_model!(m)
prepare_model!(m)

fig = Figure(resolution=(500, 800))
f1 = m.fields[1]; g1 = m.grid[1]
x1 = @view g1.x[g1.is:g1.ie]
y1 = @view g1.y[g1.js:g1.je]
s_bot1 = @view f1.s_bot[g1.is:g1.ie, g1.js:g1.je]
node1 = Node(s_bot1 .- mean(s_bot1))
h1 = heatmap(fig[1, 1], x1, y1, node1, colorrange=(-0.4, 0.4))

f2 = m.fields[2]; g2 = m.grid[2]
x2 = @view g2.x[g2.is:g2.ie]
y2 = @view g2.y[g2.js:g2.je]
s_bot2 = @view f2.s_bot[g2.is:g2.ie, g2.js:g2.je]
node2 = Node(s_bot2 .- mean(s_bot2))
h2 = heatmap(fig[2, 1], x2, y2, node2, colorrange=(-0.4, 0.4))
fig

## Run the model.
in_progress = true
while in_progress
    global in_progress = step_model!(m)
    node1[] = s_bot1 .- mean(s_bot1)
    node2[] = s_bot2 .- mean(s_bot2)
end
