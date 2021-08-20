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


# Global model data.
struct Model
    name::String
    n_domains::Int
    last_measured_time::Ref{UInt64}

    grid::Vector{Grid}
    fields::Vector{Fields}
    boundary::Vector{Boundary}
    timeloop::Vector{Timeloop}
    pressure::Vector{Pressure}
end

function Model(name, n_domains, settings)
    m = Model(name, n_domains, 0, [], [], [], [], [])
    for i in 1:n_domains
        push!(m.grid, Grid(settings[i]["grid"]))
        push!(m.fields, Fields(m.grid[i], settings[i]["fields"]))
        push!(m.boundary, Boundary(settings[i]["boundary"]))
        push!(m.timeloop, Timeloop(settings[i]["timeloop"]))
        push!(m.pressure, Pressure(m.grid[i]))
    end

    return m
end

function calc_rhs!(m::Model, i)
    set_boundary!(m.fields[i], m.grid[i], m.boundary[i])
    calc_dynamics_tend!(m.fields[i], m.grid[i])
    calc_pressure_tend!(m.fields[i], m.grid[i], m.timeloop[i], m.pressure[i])
end

function prepare_model!(m::Model)
    m.last_measured_time[] = time_ns()

    for i in 1:m.n_domains
        calc_rhs!(m, i)
    end
    check_model(m)
end

function save_domain(m::Model, i)
    f = m.fields[i]
    g = m.grid[i]
    t = m.timeloop[i]

    filename = @sprintf("%s.%02i.%08i.h5", m.name, i, round(t.time))
    h5open(filename, "w") do fid
        write(fid, "u", f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "v", f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "w", f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh])
        write(fid, "s", f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "s_bot", f.s_bot[g.is:g.ie, g.js:g.je])
        write(fid, "s_gradbot", f.s_gradbot[g.is:g.ie, g.js:g.je])
    end
end

function save_model(m::Model)
    for i in 1:m.n_domains
        save_domain(m, i)
    end
end

function load_domain!(m::Model, i)
    f = m.fields[i]
    g = m.grid[i]
    t = m.timeloop[i]

    filename = @sprintf("%s.%02i.%08i.h5", m.name, i, round(t.time))
    h5open(filename, "r") do fid
        f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "u")
        f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "v")
        f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh] = read(fid, "w")
        f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "s")
        f.s_bot[g.is:g.ie, g.js:g.je] = read(fid, "s_bot")
        f.s_gradbot[g.is:g.ie, g.js:g.je] = read(fid, "s_gradbot")
    end
end

function load_model!(m::Model)
    for i in 1:m.n_domains
        load_domain!(m, i)
    end
end

function step_model!(m::Model)
    time_next = m.timeloop[1].time + m.timeloop[1].dt

    for i in 1:m.n_domains
        while (m.timeloop[i].time < time_next)
            integrate_time!(m.fields[i], m.grid[i], m.timeloop[i])
            step_time!(m.timeloop[i])

            if (isapprox(m.timeloop[i].time % m.timeloop[i].save_time, 0.)
                && m.timeloop[i].rkstep == 1 && !isapprox(m.timeloop[i].time, m.timeloop[i].start_time))
                save_domain(m, i)
            end

            calc_rhs!(m, i)
        end
    end

    if isapprox(m.timeloop[1].time % m.timeloop[1].check_time, 0.)
        check_model(m)
    end

    return m.timeloop[1].time < m.timeloop[1].end_time
end

function check_model(m::Model)
    # Calculate the time since the last check.
    old_time = m.last_measured_time[]
    m.last_measured_time[] = time_ns()

    # First, print the model time and the wall clock since last sample.
    status_string = @sprintf("(%11.2f) Time = %8.3f",
        m.timeloop[1].time,
        (m.last_measured_time[] - old_time) * 1e-9)
    println(status_string)

    # Second, compute the divergence and CFL for each domain.
    for i in 1:m.n_domains
        status_string = @sprintf("  (%02i) Div = %6.3E, CFL = %6.3f",
            i,
            calc_divergence(m.fields[i], m.grid[i]),
            calc_cfl(m.fields[i], m.grid[i], m.timeloop[i]))
        println(status_string)
    end
end

end
