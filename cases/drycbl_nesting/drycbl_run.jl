## Loading packages.
using MicroHH
# using Profile


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings)


## Load the restart data.
load_model!(m)


## Prepare the time loop.
in_progress = prepare_model!(m)


## Start the time loop.
while in_progress
    global in_progress = step_model!(m)
    # @profile global in_progress = step_model!(m)
end

