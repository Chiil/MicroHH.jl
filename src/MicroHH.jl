module MicroHH

# Export the types.
export Model

# Export the functions.
export prepare_model!, step_model!

using LoopVectorization

include("StencilBuilder.jl")
include("Grid.jl")
include("Fields.jl")
include("Boundary.jl")
include("Timeloop.jl")
include("Dynamics.jl")
include("Pressure.jl")

struct Model
    grid::Grid
    fields::Fields
    boundary::Boundary
    timeloop::Timeloop
    pressure::Pressure
end

function Model(settings::Dict)
    grid = Grid(settings["grid"])
    fields = Fields(grid, settings["fields"])
    boundary = Boundary(settings["boundary"])
    timeloop = Timeloop(settings["timeloop"])
    pressure = Pressure(grid)

    Model(grid, fields, boundary, timeloop, pressure)
end

function calc_rhs!(model::Model)
    set_boundary!(model.fields, model.grid, model.boundary)
    calc_dynamics_tend!(model.fields, model.grid)
    calc_pressure_tend!(model.fields, model.grid, model.timeloop, model.pressure)
end

function prepare_model!(model::Model)
    calc_rhs!(model)
end

function step_model!(model::Model)
    time_next = model.timeloop.time + model.timeloop.dt
    while (model.timeloop.time < time_next)
        integrate_time!(model.fields, model.grid, model.timeloop)
        step_time!(model.timeloop)
        calc_rhs!(model)
    end

    return model.timeloop.time < model.timeloop.end_time
end

end
