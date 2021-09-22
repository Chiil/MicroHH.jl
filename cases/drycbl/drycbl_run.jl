## Loading packages.
using MicroHH
# using Profile


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings, Float32)


## Load the restart data.
load_model!(m)


## Run the model.
in_progress = prepare_model!(m)

while in_progress
    global in_progress = step_model!(m)
    # @profile global in_progress = step_model!(m)
end

