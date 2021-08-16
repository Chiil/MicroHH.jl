module MicroHH

# Export the types.
export Model

# Export the functions.
export prepare_model!, step_model!

using LoopVectorization
using Printf
using HDF5

include("StencilBuilder.jl")
include("Grid.jl")
include("Fields.jl")
include("Boundary.jl")
include("Timeloop.jl")
include("Dynamics.jl")
include("Pressure.jl")
include("Diagnostics.jl")

struct Model
    name::String
    grid::Grid
    fields::Fields
    boundary::Boundary
    timeloop::Timeloop
    pressure::Pressure
end

function Model(name, settings::Dict)
    grid = Grid(settings["grid"])
    fields = Fields(grid, settings["fields"])
    boundary = Boundary(settings["boundary"])
    timeloop = Timeloop(settings["timeloop"])
    pressure = Pressure(grid)

    Model(name, grid, fields, boundary, timeloop, pressure)
end

function calc_rhs!(model::Model)
    set_boundary!(model.fields, model.grid, model.boundary)
    calc_dynamics_tend!(model.fields, model.grid)
    calc_pressure_tend!(model.fields, model.grid, model.timeloop, model.pressure)
end

function prepare_model!(model::Model)
    save_model(model)
    calc_rhs!(model)
    println(
        "Div = ", calc_divergence(model.fields, model.grid),
        ", CFL = ", calc_cfl(model.fields, model.grid, model.timeloop))
end

function save_model(model::Model)
    f = model.fields
    g = model.grid

    filename = @sprintf("%s.%02i.%08i.h5", model.name, 1, round(model.timeloop.time))
    h5open(filename, "w") do fid
        write(fid, "u", f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "v", f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "w", f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh])
        write(fid, "s", f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
    end
end

function step_model!(model::Model)
    time_next = model.timeloop.time + model.timeloop.dt
    while (model.timeloop.time < time_next)
        integrate_time!(model.fields, model.grid, model.timeloop)
        step_time!(model.timeloop)
        if isapprox(model.timeloop.time % model.timeloop.save_time, 0.) && model.timeloop.rkstep == 1
            save_model(model)
        end
        calc_rhs!(model)
    end

    println(
        "Div = ", calc_divergence(model.fields, model.grid),
        ", CFL = ", calc_cfl(model.fields, model.grid, model.timeloop))
    return model.timeloop.time < model.timeloop.end_time
end

end
