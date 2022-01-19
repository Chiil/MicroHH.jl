module MicroHH

const do_mpi = true
const npx = 2; const npy = 2

## Export types and functions.
export Model
export prepare_model!, step_model!, save_model, load_model!


## Load packages.
if do_mpi
    using MPI
end

using LoopVectorization
using Printf
using HDF5


## Include the necessary files.
include("Parallel.jl")
include("StencilBuilder.jl")
include("Grid.jl")
include("Fields.jl")
include("Boundary.jl")
include("Timeloop.jl")
include("Dynamics.jl")
include("Pressure.jl")
include("Transposes.jl")
include("Diagnostics.jl")
include("MultiDomain.jl")


## Global model data.
struct Model{TF <: Union{Float32, Float64}}
    name::String
    n_domains::Int
    last_measured_time::Ref{UInt64}
    parallel::Parallel

    grid::Vector{Grid}
    fields::Vector{Fields}
    boundary::Vector{Boundary}
    timeloop::Vector{Timeloop}
    pressure::Vector{Pressure}
    multidomain::Vector{MultiDomain}
end


function Model(name, n_domains, settings, TF)
    parallel = Parallel(npx, npy)
    m = Model{TF}(name, n_domains, 0, parallel, [], [], [], [], [], [])
    for i in 1:n_domains
        push!(m.grid, Grid(settings[i]["grid"], TF, m.parallel))
        push!(m.fields, Fields(m.grid[i], settings[i]["fields"], TF))
        push!(m.boundary, Boundary(settings[i]["boundary"]))
        push!(m.timeloop, Timeloop(settings[i]["timeloop"]))
        push!(m.pressure, Pressure(m.grid[i], TF))
        push!(m.multidomain, MultiDomain(m.grid[i], settings[i]["multidomain"], TF))
    end

    return m
end


function calc_rhs!(m::Model, i)
    set_boundary!(m.fields[i], m.grid[i], m.boundary[i], m.parallel)
    calc_dynamics_tend!(m.fields[i], m.grid[i])
    calc_nudge_tend!(m.fields[i], m.grid[i], m.multidomain[i])
    calc_pressure_tend!(m.fields[i], m.grid[i], m.timeloop[i], m.pressure[i], m.parallel)
end


function prepare_model!(m::Model)
    # Prevent slow single precision performance due to subnormals.
    set_zero_subnormals(true)

    m.last_measured_time[] = time_ns()

    for i in 1:m.n_domains
        calc_rhs!(m, i)
    end

    check_model(m)

    return is_in_progress(m.timeloop[1])
end


function save_domain(m::Model, p::ParallelSerial)
    f = m.fields[i]
    g = m.grid[i]
    t = m.timeloop[i]
    b = m.boundary[i]

    # Update the boundary.
    set_boundary!(f, g, b, p)

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
        HDF5.h5ds_set_scale(fid["x" ], "x" )
        HDF5.h5ds_set_scale(fid["xh"], "xh")
        HDF5.h5ds_set_scale(fid["y" ], "y" )
        HDF5.h5ds_set_scale(fid["yh"], "yh")
        HDF5.h5ds_set_scale(fid["z" ], "z" )
        HDF5.h5ds_set_scale(fid["zh"], "zh")

        # Save the fields.
        write(fid, "u", f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "v", f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "w", f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh])
        write(fid, "s", f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ])
        write(fid, "s_bot", f.s_bot[g.is:g.ie, g.js:g.je])
        write(fid, "s_top", f.s_top[g.is:g.ie, g.js:g.je])
        write(fid, "s_gradbot", f.s_gradbot[g.is:g.ie, g.js:g.je])
        write(fid, "s_gradtop", f.s_gradtop[g.is:g.ie, g.js:g.je])

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

        HDF5.h5ds_attach_scale(fid["s_bot"], fid["x"], 1)
        HDF5.h5ds_attach_scale(fid["s_bot"], fid["y"], 0)

        HDF5.h5ds_attach_scale(fid["s_top"], fid["x"], 1)
        HDF5.h5ds_attach_scale(fid["s_top"], fid["y"], 0)

        HDF5.h5ds_attach_scale(fid["s_gradbot"], fid["x"], 1)
        HDF5.h5ds_attach_scale(fid["s_gradbot"], fid["y"], 0)

        HDF5.h5ds_attach_scale(fid["s_gradtop"], fid["x"], 1)
        HDF5.h5ds_attach_scale(fid["s_gradtop"], fid["y"], 0)
    end
end


function save_domain(m::Model, i, p::ParallelDistributed)
    f = m.fields[i]
    g = m.grid[i]
    t = m.timeloop[i]
    b = m.boundary[i]

    # Update the boundary.
    set_boundary!(f, g, b, p)

    # Set the range of i and j of the total grid that is stored on this task.
    is = p.id_x*g.imax + 1; ie = (p.id_x+1)*g.imax
    js = p.id_y*g.jmax + 1; je = (p.id_y+1)*g.jmax

    filename = @sprintf("%s.%02i.%08i.h5", m.name, i, round(t.time))

    # Save the grid.
    fid = h5open(filename, "w", p.commxy)
    x_id = create_dataset(fid, "x", datatype(eltype(g.x)), dataspace((g.itot,)))
    if p.id_y == 0
        x_id[is:ie] = g.x[g.is:g.ie]
    end
    MPI.Barrier(p.commx)

    xh_id = create_dataset(fid, "xh", datatype(eltype(g.xh)), dataspace((g.itot,)))
    if p.id_y == 0
        xh_id[is:ie] = g.xh[g.is:g.ie]
    end
    MPI.Barrier(p.commxy)


    y_id = create_dataset(fid, "y", datatype(eltype(g.y)), dataspace((g.jtot,)))
    if p.id_x == 0
        y_id[js:je] = g.y[g.js:g.je]
    end
    MPI.Barrier(p.commxy)

    yh_id = create_dataset(fid, "yh", datatype(eltype(g.yh)), dataspace((g.jtot,)))
    if p.id_x == 0
        yh_id[js:je] = g.yh[g.js:g.je]
    end
    MPI.Barrier(p.commxy)

    z_id = create_dataset(fid, "z", datatype(eltype(g.z)), dataspace((g.ktot,)))
    if p.id == 0
        z_id[:] = g.z[g.ks:g.ke]
    end
    MPI.Barrier(p.commxy)

    zh_id = create_dataset(fid, "zh", datatype(eltype(g.zh)), dataspace((g.ktoth,)))
    if p.id == 0
        zh_id[:] = g.zh[g.ks:g.keh]
    end
    MPI.Barrier(p.commxy)

    # Make the grid variables dimensions.
    HDF5.h5ds_set_scale(fid["x" ], "x" )
    HDF5.h5ds_set_scale(fid["xh"], "xh")
    HDF5.h5ds_set_scale(fid["y" ], "y" )
    HDF5.h5ds_set_scale(fid["yh"], "yh")
    HDF5.h5ds_set_scale(fid["z" ], "z" )
    HDF5.h5ds_set_scale(fid["zh"], "zh")

    MPI.Barrier(p.commxy)

    u_id = create_dataset(fid, "u", datatype(eltype(f.u)), dataspace((g.itot, g.jtot, g.ktot )), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    u_id[is:ie, js:je, :] = f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ]
    MPI.Barrier(p.commxy)

    v_id = create_dataset(fid, "v", datatype(eltype(f.v)), dataspace((g.itot, g.jtot, g.ktot )), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    v_id[is:ie, js:je, :] = f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ]
    MPI.Barrier(p.commxy)

    w_id = create_dataset(fid, "w", datatype(eltype(f.w)), dataspace((g.itot, g.jtot, g.ktoth)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    w_id[is:ie, js:je, :] = f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh]
    MPI.Barrier(p.commxy)

    s_id = create_dataset(fid, "s", datatype(eltype(f.s)), dataspace((g.itot, g.jtot, g.ktot )), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    s_id[is:ie, js:je, :] = f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ]
    MPI.Barrier(p.commxy)

    s_bot_id = create_dataset(fid, "s_bot", datatype(eltype(f.s_bot)), dataspace((g.itot, g.jtot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    s_bot_id[is:ie, js:je] = f.s_bot[g.is:g.ie, g.js:g.je]
    MPI.Barrier(p.commxy)

    s_top_id = create_dataset(fid, "s_top", datatype(eltype(f.s_top)), dataspace((g.itot, g.jtot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    s_top_id[is:ie, js:je] = f.s_top[g.is:g.ie, g.js:g.je]
    MPI.Barrier(p.commxy)

    s_gradbot_id = create_dataset(fid, "s_gradbot", datatype(eltype(f.s_gradbot)), dataspace((g.itot, g.jtot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    s_gradbot_id[is:ie, js:je] = f.s_gradbot[g.is:g.ie, g.js:g.je]
    MPI.Barrier(p.commxy)

    s_gradtop_id = create_dataset(fid, "s_gradtop", datatype(eltype(f.s_gradtop)), dataspace((g.itot, g.jtot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    s_gradtop_id[is:ie, js:je] = f.s_gradtop[g.is:g.ie, g.js:g.je]
    MPI.Barrier(p.commxy)

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

    HDF5.h5ds_attach_scale(fid["s_bot"], fid["x"], 1)
    HDF5.h5ds_attach_scale(fid["s_bot"], fid["y"], 0)

    HDF5.h5ds_attach_scale(fid["s_top"], fid["x"], 1)
    HDF5.h5ds_attach_scale(fid["s_top"], fid["y"], 0)

    HDF5.h5ds_attach_scale(fid["s_gradbot"], fid["x"], 1)
    HDF5.h5ds_attach_scale(fid["s_gradbot"], fid["y"], 0)

    HDF5.h5ds_attach_scale(fid["s_gradtop"], fid["x"], 1)
    HDF5.h5ds_attach_scale(fid["s_gradtop"], fid["y"], 0)

    close(fid)
end


function save_model(m::Model)
    for i in 1:m.n_domains
        save_domain(m, i, m.parallel)
    end
end


function load_domain!(m::Model, i)
    f = m.fields[i]
    g = m.grid[i]
    t = m.timeloop[i]
    b = m.boundary[i]
    p = m.parallel

    filename = @sprintf("%s.%04i.%02i.%08i.h5", m.name, m.parallel.id, i, round(t.time))
    h5open(filename, "r") do fid
        f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "u")
        f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "v")
        f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh] = read(fid, "w")
        f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "s")
        f.s_bot[g.is:g.ie, g.js:g.je] = read(fid, "s_bot")
        f.s_top[g.is:g.ie, g.js:g.je] = read(fid, "s_top")
        f.s_gradbot[g.is:g.ie, g.js:g.je] = read(fid, "s_gradbot")
        f.s_gradtop[g.is:g.ie, g.js:g.je] = read(fid, "s_gradtop")
    end

    set_boundary!(f, g, b, p)
end


function load_model!(m::Model)
    for i in 1:m.n_domains
        load_domain!(m, i)
    end
end


function step_model!(m::Model)
    # Prevent slow single precision performance due to subnormals.
    set_zero_subnormals(true)

    itime_next = m.timeloop[1].itime + m.timeloop[1].idt

    for i in 1:m.n_domains
        # Compute the nudging fields, which remain constant throughout the step.
        if i > 1
            calc_nudge_fields!(m.multidomain[i], m.fields[i], m.fields[i-1], m.grid[i], m.grid[i-1])
        end

        while (m.timeloop[i].itime < itime_next)
            integrate_time!(m.fields[i], m.grid[i], m.timeloop[i])
            step_time!(m.timeloop[i])

            if is_save_time(m.timeloop[i])
                save_domain(m, i, m.parallel)
            end

            calc_rhs!(m, i)
        end
    end

    if is_check_time(m.timeloop[1])
        check_model(m)
    end

    return is_in_progress(m.timeloop[1])
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
        status_string = @sprintf("  (%02i) Div = %6.3E, CFL = %6.3f, subnormals = %i",
            i,
            calc_divergence(m.fields[i], m.grid[i]),
            calc_cfl(m.fields[i], m.grid[i], m.timeloop[i]),
            sum(issubnormal.(m.fields[i].u)) + sum(issubnormal.(m.fields[i].v)) + sum(issubnormal.(m.fields[i].w)))
        println(status_string)
    end
end

end
