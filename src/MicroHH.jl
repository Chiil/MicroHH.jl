module MicroHH

export Model, Grid, Fields, Dynamics

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

end
