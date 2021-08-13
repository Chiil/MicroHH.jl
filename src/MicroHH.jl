module MicroHH

# Export the types.
export Model

# Export the functions.
export step_model!

using LoopVectorization

include("StencilBuilder.jl")
include("Grid.jl")
include("Fields.jl")
include("Boundary.jl")
include("Timeloop.jl")
include("Dynamics.jl")

struct Model
    grid::Grid
    fields::Fields
    boundary::Boundary
    timeloop::Timeloop
end

function Model(settings::Dict)
    grid = Grid(settings["grid"])
    fields = Fields(grid, settings["fields"])
    boundary = Boundary(settings["boundary"])
    timeloop = Timeloop(settings["timeloop"])

    Model(grid, fields, boundary, timeloop)
end

function step_model!(model::Model)
    set_boundary!(model.fields, model.grid, model.boundary)
    calc_dynamics!(model.fields, model.grid)
    return step_time!(model.timeloop)
end

end
