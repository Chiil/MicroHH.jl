module MicroHH

## Exports
# Export the types.
export Model

# Export the functions.
export prepare_model!, step_model!, save_model, load_model!

using LoopVectorization
using Printf
using HDF5


## Prevent slow single precision performance due to subnormals.
set_zero_subnormals(true)


## Include the necessary files.
include("StencilBuilder.jl")
include("Grid.jl")
include("Fields.jl")
include("Boundary.jl")
include("Timeloop.jl")
include("Dynamics.jl")
include("Pressure.jl")
include("Diagnostics.jl")
include("MultiDomain.jl")


## Global model data.
struct Model{TF <: Union{Float32, Float64}}
    name::String
    n_domains::Int
    last_measured_time::Ref{UInt64}

    grid::Vector{Grid}
    fields::Vector{Fields}
    boundary::Vector{Boundary}
    timeloop::Vector{Timeloop}
    pressure::Vector{Pressure}
    multidomain::Vector{MultiDomain}
end


function Model(name, n_domains, settings, TF)
    m = Model{TF}(name, n_domains, 0, [], [], [], [], [], [])
    for i in 1:n_domains
        push!(m.grid, Grid(settings[i]["grid"], TF))
        push!(m.fields, Fields(m.grid[i], settings[i]["fields"], TF))
        push!(m.boundary, Boundary(settings[i]["boundary"]))
        push!(m.timeloop, Timeloop(settings[i]["timeloop"]))
        push!(m.pressure, Pressure(m.grid[i], TF))
        push!(m.multidomain, MultiDomain(m.grid[i], settings[i]["multidomain"], TF))
    end

    return m
end


function calc_rhs!(m::Model, i)
    set_boundary!(m.fields[i], m.grid[i], m.boundary[i])
    calc_dynamics_tend!(m.fields[i], m.grid[i])
    calc_nudge_tend!(m.fields[i], m.grid[i], m.multidomain[i])
    calc_pressure_tend!(m.fields[i], m.grid[i], m.timeloop[i], m.pressure[i])
end


function prepare_model!(m::Model)
    m.last_measured_time[] = time_ns()

    for i in 1:m.n_domains
        calc_rhs!(m, i)
    end

    check_model(m)

    return m.timeloop[1].itime < m.timeloop[1].iend_time
end


function save_domain(m::Model, i)
    f = m.fields[i]
    g = m.grid[i]
    t = m.timeloop[i]

    filename = @sprintf("%s.%02i.%08i.h5", m.name, i, round(t.time))
    h5open(filename, "w") do fid
        # Save the grid.
        write(fid, "x" , g.x[g.is:g.ie])
        write(fid, "xh", g.xh[g.is:g.ie])
        write(fid, "y" , g.y[g.js:g.je])
        write(fid, "yh", g.yh[g.js:g.je])
        write(fid, "z" , g.z[g.ks:g.ke])
        write(fid, "zh", g.zh[g.ks:g.keh])

        # Make the grid variables dimensions.
        HDF5.h5ds_set_scale(fid["x"], "x")
        HDF5.h5ds_set_scale(fid["xh"], "xh")
        HDF5.h5ds_set_scale(fid["y"], "y")
        HDF5.h5ds_set_scale(fid["yh"], "yh")
        HDF5.h5ds_set_scale(fid["z"], "z")
        HDF5.h5ds_set_scale(fid["zh"], "zh")

        # Save the fields.
        write(fid, "u", f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "v", f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "w", f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh])
        write(fid, "s", f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "s_bot", f.s_bot[g.is:g.ie, g.js:g.je])
        write(fid, "s_gradbot", f.s_gradbot[g.is:g.ie, g.js:g.je])

        # Attach the dimensions. Note the c-indexing.
        HDF5.h5ds_attach_scale(fid["u"], fid["xh"], 2)
        HDF5.h5ds_attach_scale(fid["u"], fid["y"], 1)
        HDF5.h5ds_attach_scale(fid["u"], fid["z"], 0)

        HDF5.h5ds_attach_scale(fid["v"], fid["x"], 2)
        HDF5.h5ds_attach_scale(fid["v"], fid["yh"], 1)
        HDF5.h5ds_attach_scale(fid["v"], fid["z"], 0)

        HDF5.h5ds_attach_scale(fid["w"], fid["x"], 2)
        HDF5.h5ds_attach_scale(fid["w"], fid["y"], 1)
        HDF5.h5ds_attach_scale(fid["w"], fid["zh"], 0)

        HDF5.h5ds_attach_scale(fid["s"], fid["x"], 2)
        HDF5.h5ds_attach_scale(fid["s"], fid["y"], 1)
        HDF5.h5ds_attach_scale(fid["s"], fid["z"], 0)
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
    itime_next = m.timeloop[1].itime + m.timeloop[1].idt

    for i in 1:m.n_domains
        # Compute the nudging fields, which remain constant throughout the step.
        if i > 1
            calc_nudge_fields!(m.multidomain[i], m.fields[i], m.fields[i-1], m.grid[i], m.grid[i-1])
        end

        while (m.timeloop[i].itime < itime_next)
            integrate_time!(m.fields[i], m.grid[i], m.timeloop[i])
            step_time!(m.timeloop[i])

            if (m.timeloop[i].itime % m.timeloop[i].isave_time == 0
                && m.timeloop[i].rkstep == 1 && m.timeloop[i].itime != m.timeloop[i].istart_time)
                save_domain(m, i)
            end

            calc_rhs!(m, i)
        end
    end

    if m.timeloop[1].itime % m.timeloop[1].icheck_time == 0
        check_model(m)
    end

    return m.timeloop[1].itime < m.timeloop[1].iend_time
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
