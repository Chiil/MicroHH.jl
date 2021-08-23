## Loading packages.
using MicroHH
using Tullio
using Statistics
using GLMakie


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings, Float32)
load_model!(m)
prepare_model!(m)
f = m.fields[1]; g = m.grid[1]
x = @view g.x[g.is:g.ie]
y = @view g.y[g.js:g.je]
s_bot = @view f.s_bot[g.is:g.ie, g.js:g.je]
node = Node(s_bot .- mean(s_bot))
heatmap(x, y, node, colorrange=(-0.4, 0.4))


## Run the model.
in_progress = true
while in_progress
    global in_progress = step_model!(m)
    node[] = s_bot .- mean(s_bot)
end
