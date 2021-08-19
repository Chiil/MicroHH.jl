module MicroHH

# Export the types.
export Model

# Export the functions.
export prepare_model!, step_model!, save_model, load_model!

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

    last_measured_time::Ref{UInt64}
end

function Model(name, settings::Dict)
    grid = Grid(settings["grid"])
    fields = Fields(grid, settings["fields"])
    boundary = Boundary(settings["boundary"])
    timeloop = Timeloop(settings["timeloop"])
    pressure = Pressure(grid)

    Model(name, grid, fields, boundary, timeloop, pressure, 0)
end

function calc_rhs!(m::Model)
    set_boundary!(m.fields, m.grid, m.boundary)
    calc_dynamics_tend!(m.fields, m.grid)
    calc_pressure_tend!(m.fields, m.grid, m.timeloop, m.pressure)
end

function prepare_model!(m::Model)
    calc_rhs!(m)
    m.last_measured_time[] = time_ns()
    check_model(m)
end

function save_model(m::Model)
    f = m.fields
    g = m.grid

    filename = @sprintf("%s.%02i.%08i.h5", m.name, 1, round(m.timeloop.time))
    h5open(filename, "w") do fid
        write(fid, "u", f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "v", f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "w", f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh])
        write(fid, "s", f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "s_bot", f.s_bot[g.is:g.ie, g.js:g.je])
        write(fid, "s_gradbot", f.s_gradbot[g.is:g.ie, g.js:g.je])
    end
end

function load_model!(m::Model)
    f = m.fields
    g = m.grid

    filename = @sprintf("%s.%02i.%08i.h5", m.name, 1, round(m.timeloop.time))
    h5open(filename, "r") do fid
        f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "u")
        f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "v")
        f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh] = read(fid, "w")
        f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "s")
        f.s_bot[g.is:g.ie, g.js:g.je] = read(fid, "s_bot")
        f.s_gradbot[g.is:g.ie, g.js:g.je] = read(fid, "s_gradbot")
    end
end

function step_model!(m::Model)
    time_next = m.timeloop.time + m.timeloop.dt

    while (m.timeloop.time < time_next)
        integrate_time!(m.fields, m.grid, m.timeloop)
        step_time!(m.timeloop)

        if (isapprox(m.timeloop.time % m.timeloop.save_time, 0.)
            && m.timeloop.rkstep == 1 && !isapprox(m.timeloop.time, m.timeloop.start_time))
            save_model(m)
        end

        calc_rhs!(m)
    end

    if isapprox(m.timeloop.time % m.timeloop.check_time, 0.)
        check_model(m)
    end

    return m.timeloop.time < m.timeloop.end_time
end

function check_model(m::Model)
    old_time = m.last_measured_time[]
    m.last_measured_time[] = time_ns()
    status_string = @sprintf("(%11.2f) Time/iter = %8.3f, Div = %6.3E, CFL = %6.3f",
        m.timeloop.time,
        (m.last_measured_time[] - old_time) * 1e-9 * m.timeloop.dt / m.timeloop.check_time,
        calc_divergence(m.fields, m.grid),
        calc_cfl(m.fields, m.grid, m.timeloop))

    println(status_string)
end

end
