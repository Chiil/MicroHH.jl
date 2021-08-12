module MicroHH

# Export the types.
export Model

# Export the functions.
export step_model

include("Grid.jl")
include("Fields.jl")
include("Dynamics.jl")

struct Model
    grid::Grid
    fields::Fields
end

function Model(settings::Dict)
    grid = Grid(settings["grid"])
    fields = Fields(grid, settings["fields"])

    Model(grid, fields)
end

function step_model(model::Model)
    calc_dynamics!(model.fields, model.grid)
    return true
end

end
