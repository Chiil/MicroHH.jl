module MicroHH


## Export types and functions.
export Model
export prepare_model!, step_model!, save_model, load_model!, output_timer_model!


## Load packages.
const use_mpi = Ref{Bool}(false)
const npx = Ref{Int}(1); const npy = Ref{Int}(1)

using LoopVectorization
using Tullio
using Printf
using HDF5
using ArgParse
using TimerOutputs


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
    to::TimerOutput

    grid::Vector{Grid}
    fields::Vector{Fields}
    boundary::Vector{Boundary}
    timeloop::Vector{Timeloop}
    pressure::Vector{Pressure}
    multidomain::Vector{MultiDomain}
end


function Model(name, n_domains, settings, TF)
    parallel = Parallel(npx[], npy[])
    to = TimerOutput()
    disable_timer!(to)

    m = Model{TF}(name, n_domains, 0, parallel, to, [], [], [], [], [], [])
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
    @timeit m.to "set_boundary"       set_boundary!(m.fields[i], m.grid[i], m.boundary[i], m.parallel)
    @timeit m.to "calc_dynamics_tend" calc_dynamics_tend!(m.fields[i], m.grid[i], m.to)
    @timeit m.to "calc_nudge_tend"    calc_nudge_tend!(m.fields[i], m.grid[i], m.multidomain[i])
    @timeit m.to "calc_pressure_tend" calc_pressure_tend!(m.fields[i], m.grid[i], m.timeloop[i], m.pressure[i], m.parallel, m.to)
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


function save_domain(m::Model, i, p::ParallelSerial)
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
        write(fid, "s_ref", f.s_ref[g.ks:g.ke])

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

        HDF5.h5ds_attach_scale(fid["s_ref"], fid["z"], 0)
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
    close(x_id)

    xh_id = create_dataset(fid, "xh", datatype(eltype(g.xh)), dataspace((g.itot,)))
    if p.id_y == 0
        xh_id[is:ie] = g.xh[g.is:g.ie]
    end
    close(xh_id)

    y_id = create_dataset(fid, "y", datatype(eltype(g.y)), dataspace((g.jtot,)))
    if p.id_x == 0
        y_id[js:je] = g.y[g.js:g.je]
    end
    close(y_id)

    yh_id = create_dataset(fid, "yh", datatype(eltype(g.yh)), dataspace((g.jtot,)))
    if p.id_x == 0
        yh_id[js:je] = g.yh[g.js:g.je]
    end
    close(yh_id)

    z_id = create_dataset(fid, "z", datatype(eltype(g.z)), dataspace((g.ktot,)))
    if p.id == 0
        z_id[:] = g.z[g.ks:g.ke]
    end
    close(z_id)

    zh_id = create_dataset(fid, "zh", datatype(eltype(g.zh)), dataspace((g.ktoth,)))
    if p.id == 0
        zh_id[:] = g.zh[g.ks:g.keh]
    end
    close(zh_id)

    # Make the grid variables dimensions.
    HDF5.h5ds_set_scale(fid["x" ], "x" )
    HDF5.h5ds_set_scale(fid["xh"], "xh")
    HDF5.h5ds_set_scale(fid["y" ], "y" )
    HDF5.h5ds_set_scale(fid["yh"], "yh")
    HDF5.h5ds_set_scale(fid["z" ], "z" )
    HDF5.h5ds_set_scale(fid["zh"], "zh")

    MPI.Barrier(p.commxy)

    items_to_save = [
        ("u", f.u, g.ktot ),
        ("v", f.v, g.ktot ),
        ("w", f.w, g.ktoth),
        ("s", f.s, g.ktot )]

    map(items_to_save) do item
        name, a, ktot = item
        a_nogc = a[g.is:g.ie, g.js:g.je, g.ks:g.ks+ktot-1]
        aid = create_dataset(fid, name, datatype(a_nogc), dataspace((g.itot, g.jtot, ktot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
        aid[is:ie, js:je, :] = a_nogc[:, :, :]
        close(aid)
    end

    items_to_save = [
        ("s_bot", f.s_bot),
        ("s_top", f.s_top),
        ("s_gradbot", f.s_gradbot),
        ("s_gradtop", f.s_gradtop)]

    map(items_to_save) do item
        name, a = item
        a_nogc = a[g.is:g.ie, g.js:g.je]
        aid = create_dataset(fid, name, datatype(a_nogc), dataspace((g.itot, g.jtot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
        aid[is:ie, js:je] = a_nogc[:, :]
        close(aid)
    end

    s_ref_id = create_dataset(fid, "s_ref", datatype(eltype(f.s_ref)), dataspace((g.ktot,)))
    if p.id == 0
        s_ref_id[:] = f.s_ref[g.ks:g.ke]
    end
    close(s_ref_id)

    # Attach the dimensions. Note the c-indexing.
    vars_3d = [
        ("u", "xh", "y" , "z" ),
        ("v", "x" , "yh", "z" ),
        ("w", "x" , "y" , "zh"),
        ("s", "x" , "y" , "z" )]

    map(vars_3d) do var
        name::String, x::String, y::String, z::String = var
        HDF5.h5ds_attach_scale(fid[name], fid[x], 2)
        HDF5.h5ds_attach_scale(fid[name], fid[y], 1)
        HDF5.h5ds_attach_scale(fid[name], fid[z], 0)
    end

    vars_2d = [
        ("s_bot", "x", "y"),
        ("s_top", "x", "y"),
        ("s_gradbot", "x", "y"),
        ("s_gradtop", "x", "y")]

    map(vars_2d) do var
        name::String, x::String, y::String = var
        HDF5.h5ds_attach_scale(fid[name], fid[x], 1)
        HDF5.h5ds_attach_scale(fid[name], fid[y], 0)
    end

    HDF5.h5ds_attach_scale(fid["s_ref"], fid["z"], 0)

    # Close
    close(fid)
end


function save_model(m::Model)
    for i in 1:m.n_domains
        save_domain(m, i, m.parallel)
    end
end


function load_domain!(m::Model, i, p::ParallelSerial)
    f = m.fields[i]
    g = m.grid[i]
    t = m.timeloop[i]
    b = m.boundary[i]
    p = m.parallel

    filename = @sprintf("%s.%02i.%08i.h5", m.name, i, round(t.time))
    h5open(filename, "r") do fid
        f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "u")
        f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "v")
        f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh] = read(fid, "w")
        f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = read(fid, "s")
        f.s_bot[g.is:g.ie, g.js:g.je] = read(fid, "s_bot")
        f.s_top[g.is:g.ie, g.js:g.je] = read(fid, "s_top")
        f.s_gradbot[g.is:g.ie, g.js:g.je] = read(fid, "s_gradbot")
        f.s_gradtop[g.is:g.ie, g.js:g.je] = read(fid, "s_gradtop")
        f.s_ref[g.ks:g.ke] = read(fid, "s_ref")
    end

    # Set the boundary values and ghost cells.
    set_boundary!(f, g, b, p)

    # Extrapolate the reference profile to the ghost cells.
    extrap = LinearInterpolation(g.z[g.ks:g.ke], f.s_ref[g.ks:g.ke], extrapolation_bc=Line())
    f.s_ref[1:g.kgc] = extrap(g.z[1:g.kgc])
    f.s_ref[g.ke+1:end] = extrap(g.z[g.ke+1:end])
end


function load_domain!(m::Model, i, p::ParallelDistributed)
    f = m.fields[i]
    g = m.grid[i]
    t = m.timeloop[i]
    b = m.boundary[i]
    p = m.parallel

    # Set the range of i and j of the total grid that is stored on this task.
    is = p.id_x*g.imax + 1; ie = (p.id_x+1)*g.imax
    js = p.id_y*g.jmax + 1; je = (p.id_y+1)*g.jmax

    filename = @sprintf("%s.%02i.%08i.h5", m.name, i, round(t.time))
    fid = h5open(filename, "r", p.commxy)

    dapl = create_property(HDF5.H5P_DATASET_ACCESS)
    dxpl = create_property(HDF5.H5P_DATASET_XFER)
    dxpl[:dxpl_mpio] = HDF5.H5FD_MPIO_COLLECTIVE

    u_id = open_dataset(fid, "u", dapl, dxpl)
    f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = u_id[is:ie, js:je, :]
    close(u_id)

    v_id = open_dataset(fid, "v", dapl, dxpl)
    f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = v_id[is:ie, js:je, :]
    close(v_id)

    w_id = open_dataset(fid, "w", dapl, dxpl)
    f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh] = w_id[is:ie, js:je, :]
    close(w_id)

    s_id = open_dataset(fid, "s", dapl, dxpl)
    f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ] = s_id[is:ie, js:je, :]
    close(s_id)

    s_bot_id = open_dataset(fid, "s_bot", dapl, dxpl)
    f.s_bot[g.is:g.ie, g.js:g.je] = s_bot_id[is:ie, js:je]
    close(s_bot_id)

    s_top_id = open_dataset(fid, "s_top", dapl, dxpl)
    f.s_top[g.is:g.ie, g.js:g.je] = s_top_id[is:ie, js:je]
    close(s_top_id)

    s_gradbot_id = open_dataset(fid, "s_gradbot", dapl, dxpl)
    f.s_gradbot[g.is:g.ie, g.js:g.je] = s_gradbot_id[is:ie, js:je]
    close(s_gradbot_id)

    s_gradtop_id = open_dataset(fid, "s_gradtop", dapl, dxpl)
    f.s_gradtop[g.is:g.ie, g.js:g.je] = s_gradtop_id[is:ie, js:je]
    close(s_gradtop_id)

    s_ref_id = open_dataset(fid, "s_ref", dapl, dxpl)
    f.s_ref[g.ks:g.ke] = s_ref_id[:]
    close(s_ref_id)

    close(dapl)
    close(dxpl)
    close(fid)

    # Set the boundary values and ghost cells.
    set_boundary!(f, g, b, p)

    # Extrapolate the reference profile to the ghost cells.
    extrap = LinearInterpolation(g.z[g.ks:g.ke], f.s_ref[g.ks:g.ke], extrapolation_bc=Line())
    f.s_ref[1:g.kgc] = extrap(g.z[1:g.kgc])
    f.s_ref[g.ke+1:end] = extrap(g.z[g.ke+1:end])
end


function load_model!(m::Model)
    for i in 1:m.n_domains
        load_domain!(m, i, m.parallel)
    end
end


function step_model!(m::Model)
    # Prevent slow single precision performance due to subnormals.
    set_zero_subnormals(true)

    itime_next = m.timeloop[1].itime + m.timeloop[1].idt

    @timeit m.to "domains" begin
        for i in 1:m.n_domains
            # Compute the nudging fields, which remain constant throughout the step.
            if i > 1
                @timeit m.to "calc_nudge_fields" calc_nudge_fields!(m.multidomain[i], m.fields[i], m.fields[i-1], m.grid[i], m.grid[i-1])
            end

            while (m.timeloop[i].itime < itime_next)
                @timeit m.to "integrate_time" integrate_time!(m.fields[i], m.grid[i], m.timeloop[i])
                @timeit m.to "step_time" step_time!(m.timeloop[i])

                if is_save_time(m.timeloop[i])
                    @timeit m.to "save_domain" save_domain(m, i, m.parallel)
                end

                @timeit m.to "calc_rhs" calc_rhs!(m, i)
            end
        end
    end

    if is_check_time(m.timeloop[1])
        @timeit m.to "check_model" check_model(m)
    end

    # Enable the timer after the first round to avoid measuring compilation.
    enable_timer!(m.to)

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
    if m.parallel.id == 0
        @info "$status_string"
    end

    # Second, compute the divergence and CFL for each domain.
    for i in 1:m.n_domains
        status_string = @sprintf("  (%02i) Div = %6.3E, CFL = %6.3f, subnormals = %i",
            i,
            calc_divergence(m.fields[i], m.grid[i]),
            calc_cfl(m.fields[i], m.grid[i], m.timeloop[i]),
            sum(issubnormal.(m.fields[i].u)) + sum(issubnormal.(m.fields[i].v)) + sum(issubnormal.(m.fields[i].w)))

        if m.parallel.id == 0
            @info "$status_string"
        end
    end
end


function output_timer_model!(m::Model)
    if m.parallel.id == 0
        @info ""
        @info "Timer output:"
        show(m.to)
    end
end


function __init__()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--use-mpi"
            help = "Enable MPI"
            action = :store_true
        "--npx"
            help = "MPI pencils in x-direction"
            arg_type = Int
            default = 1
        "--npy"
            help = "MPI pencils in y-direction"
            arg_type = Int
            default = 1
    end

    parsed_args = parse_args(s)

    use_mpi[] = parsed_args["use-mpi"]
    npx[] = parsed_args["npx"]
    npy[] = parsed_args["npy"]

    # Only load MPI if parallel run is required.
    if use_mpi[]
        @eval using MPI
    end
end


end
