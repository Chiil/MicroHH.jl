## Loading packages.
using MicroHH


## Loading settings.
include("lid_driven_cavity_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("lid_driven_cavity", n_domains, settings, float_type)


## Load the restart data.
load_model!(m)


## Run the model.
in_progress = prepare_model!(m)

while in_progress
    global in_progress = step_model!(m)
end

output_timer_model!(m)
