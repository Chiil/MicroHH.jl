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


## Run the model.
prepare_model!(m)

in_progress = true
while in_progress
    global in_progress = step_model!(m)
    # @profile global in_progress = step_model!(m)
end

