module MicroHH

# Export the types.
export Model

# Export the functions.
export step_model

using LoopVectorization

include("StencilBuilder.jl")
include("Grid.jl")
include("Fields.jl")
include("Dynamics.jl")
include("Timeloop.jl")

struct Model
    grid::Grid
    fields::Fields
    timeloop::Timeloop
end

function Model(settings::Dict)
    grid = Grid(settings["grid"])
    fields = Fields(grid, settings["fields"])
    timeloop = Timeloop(settings["timeloop"])

    Model(grid, fields, timeloop)
end

function step_model(model::Model)
    println(model.timeloop.time)
    calc_dynamics!(model.fields, model.grid)
    return step_time!(model.timeloop)
end

end
