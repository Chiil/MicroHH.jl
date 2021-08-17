## Loading packages.
using MicroHH
using Tullio
using Statistics
# using Profile


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
model = Model("drycbl", settings)


## Load the restart data.
load_model!(model)


## Run the model.
prepare_model!(model)

in_progress = true
while in_progress
    duration = @timed global in_progress = step_model!(model)
    println(duration.time, " ", model.timeloop.time)
    # @profile global in_progress = step_model!(model)
end

