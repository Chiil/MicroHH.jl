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
heatmap(m.domains[1].fields.s[:, :, 2])


## Run the model.
in_progress = true
while in_progress
    global in_progress = step_model!(m)
    heatmap!(m.domains[1].fields.s[:, :, 2])
end
