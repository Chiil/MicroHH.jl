## Loading packages.
using MicroHH
using Tullio
using Statistics
using GLMakie


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings)
load_model!(m)
prepare_model!(m)
f = m.fields[1]; g = m.grid[1]
s = @view f.s[g.is:g.ie, g.js:g.je, 2]
node = Node(s .- mean(s))
heatmap(node)


## Run the model.
in_progress = true
while in_progress
    global in_progress = step_model!(m)
    node[] = s .- mean(s)
end
