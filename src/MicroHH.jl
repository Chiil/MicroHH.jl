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

# Domain data.
struct Domain
    grid::Grid
    fields::Fields
    boundary::Boundary
    timeloop::Timeloop
    pressure::Pressure
end

function Domain(settings::Dict)
    grid = Grid(settings["grid"])
    fields = Fields(grid, settings["fields"])
    boundary = Boundary(settings["boundary"])
    timeloop = Timeloop(settings["timeloop"])
    pressure = Pressure(grid)

    Domain(grid, fields, boundary, timeloop, pressure)
end

# Global model data.
struct Model
    name::String
    n_domains
    domains::Vector{Domain}
    last_measured_time::Ref{UInt64}
end

function Model(name, n_domains, settings)
    m = Model(name, n_domains, [], 0)
    for i in 1:n_domains
        push!(m.domains, Domain(settings[i]))
    end

    return m
end

function calc_rhs!(d::Domain)
    set_boundary!(d.fields, d.grid, d.boundary)
    calc_dynamics_tend!(d.fields, d.grid)
    calc_pressure_tend!(d.fields, d.grid, d.timeloop, d.pressure)
end

function prepare_model!(m::Model)
    m.last_measured_time[] = time_ns()

    for i in 1:m.n_domains
        d = m.domains[i]
        calc_rhs!(d)
    end
    check_model(m)
end

function save_domain(d::Domain, name, i_domain)
    f = d.fields
    g = d.grid
    t = d.timeloop

    filename = @sprintf("%s.%02i.%08i.h5", name, i_domain, round(t.time))
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
        save_domain(m.domains[i], m.name, i)
    end
end

function load_domain!(d::Domain, name, i_domain)
    f = d.fields
    g = d.grid
    t = d.timeloop

    filename = @sprintf("%s.%02i.%08i.h5", name, i_domain, round(t.time))
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
        load_domain!(m.domains[i], m.name, i)
    end
end

function step_model!(m::Model)
    time_next = m.domains[1].timeloop.time + m.domains[1].timeloop.dt

    for i in 1:m.n_domains
        d = m.domains[i]

        while (d.timeloop.time < time_next)
            integrate_time!(d.fields, d.grid, d.timeloop)
            step_time!(d.timeloop)

            if (isapprox(d.timeloop.time % d.timeloop.save_time, 0.)
                && d.timeloop.rkstep == 1 && !isapprox(d.timeloop.time, d.timeloop.start_time))
                save_domain(d, m.name, i)
            end

            calc_rhs!(d)
        end
    end

    if isapprox(m.domains[1].timeloop.time % m.domains[1].timeloop.check_time, 0.)
        check_model(m)
    end

    return m.domains[1].timeloop.time < m.domains[1].timeloop.end_time
end

function check_model(m::Model)
    old_time = m.last_measured_time[]
    m.last_measured_time[] = time_ns()

    # First, print the model time and the wall clock since last sample.
    status_string = @sprintf("(%11.2f) Time = %8.3f",
        m.domains[1].timeloop.time,
        (m.last_measured_time[] - old_time) * 1e-9)
    println(status_string)

    # Second, compute the divergence and CFL for each domain.
    for i in 1:m.n_domains
        d = m.domains[i]
        status_string = @sprintf("  (%02i) Div = %6.3E, CFL = %6.3f",
            i,
            calc_divergence(d.fields, d.grid),
            calc_cfl(d.fields, d.grid, d.timeloop))
        println(status_string)
    end
end

end
