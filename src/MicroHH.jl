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
include("Stats.jl")


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
        push!(m.grid, Grid(settings[i]["grid"], m.parallel, TF))
        push!(m.fields, Fields(m.grid[i], settings[i]["fields"], TF))
        push!(m.boundary, Boundary(m.grid[i], m.parallel, settings[i]["boundary"], TF))
        push!(m.timeloop, Timeloop(settings[i]["timeloop"]))
        push!(m.pressure, Pressure(m.grid[i], m.parallel, TF))
        push!(m.multidomain, MultiDomain(m.grid[i], settings[i]["multidomain"], TF))
    end

    return m
end


function calc_rhs!(m::Model, i)
    @timeit m.to "set_boundary"       set_boundary!(m.fields[i], m.grid[i], m.boundary[i], m.parallel)

    # CvH TMP
    if i > 1
        g = m.grid[i]; gsrc = m.grid[i-1]
        f = m.fields[i]; fsrc = m.fields[i-1]

        # HARDCODE THE NESTING, FIX THIS ONCE IT WORKS
        isc = gsrc.is + gsrc.itot ÷ 4; iec = gsrc.is + 3*(gsrc.itot÷4)
        jsc = gsrc.js; jec = gsrc.je
        ksc = gsrc.ks; kec = gsrc.ke

        # Set the walls to the parent domain.
        u_west = @view fsrc.u[isc, jsc-1:jec+1, ksc-1:kec+1]
        u_east = @view fsrc.u[iec, jsc-1:jec+1, ksc-1:kec+1]

        interp = interpolate((gsrc.y, gsrc.z), u_west, (Gridded(Linear()), Gridded(Linear())))
        f.u[g.is, :, :] .= interp(g.y, g.z)
        interp = interpolate((gsrc.y, gsrc.z), u_east, (Gridded(Linear()), Gridded(Linear())))
        f.u[g.ie+1, :, :] .= interp(g.y, g.z)

        v_west_p = 0.5 .* (  fsrc.v[isc-1, jsc-1:jec+1, ksc-1:kec+1]
                          .+ fsrc.v[isc  , jsc-1:jec+1, ksc-1:kec+1] )
        v_east_p = 0.5 .* (  fsrc.v[iec  , jsc-1:jec+1, ksc-1:kec+1]
                          .+ fsrc.v[iec+1, jsc-1:jec+1, ksc-1:kec+1] )

        interp = interpolate((gsrc.y, gsrc.z), v_west_p, (Gridded(Linear()), Gridded(Linear())))
        v_west = interp(g.y, g.z)
        interp = interpolate((gsrc.y, gsrc.z), v_east_p, (Gridded(Linear()), Gridded(Linear())))
        v_east = interp(g.y, g.z)

        @. f.v[g.is-1, :, :] = 2 * v_west - f.v[g.is, :, :]
        @. f.v[g.ie+1, :, :] = 2 * v_east - f.v[g.ie, :, :]

        w_west_p = 0.5 .* (  fsrc.w[isc-1, jsc-1:jec+1, ksc-1:kec+1]
                          .+ fsrc.w[isc  , jsc-1:jec+1, ksc-1:kec+1] )
        w_east_p = 0.5 .* (  fsrc.w[iec  , jsc-1:jec+1, ksc-1:kec+1]
                          .+ fsrc.w[iec+1, jsc-1:jec+1, ksc-1:kec+1] )

        interp = interpolate((gsrc.y, gsrc.zh), w_west_p, (Gridded(Linear()), Gridded(Linear())))
        w_west = interp(g.y, g.zh)
        interp = interpolate((gsrc.y, gsrc.zh), w_east_p, (Gridded(Linear()), Gridded(Linear())))
        w_east = interp(g.y, g.zh)

        @. f.w[g.is-1, :, :] = 2 * w_west - f.w[g.is, :, :]
        @. f.w[g.ie+1, :, :] = 2 * w_east - f.w[g.ie, :, :]

        s_west_p = 0.5 .* (  fsrc.s[isc-1, jsc-1:jec+1, ksc-1:kec+1]
                          .+ fsrc.s[isc  , jsc-1:jec+1, ksc-1:kec+1] )
        s_east_p = 0.5 .* (  fsrc.s[iec  , jsc-1:jec+1, ksc-1:kec+1]
                          .+ fsrc.s[iec+1, jsc-1:jec+1, ksc-1:kec+1] )

        interp = interpolate((gsrc.y, gsrc.z), s_west_p, (Gridded(Linear()), Gridded(Linear())))
        s_west = interp(g.y, g.z)
        interp = interpolate((gsrc.y, gsrc.z), s_east_p, (Gridded(Linear()), Gridded(Linear())))
        s_east = interp(g.y, g.z)

        @. f.s[g.is-1, :, :] = 2 * s_west - f.s[g.is, :, :]
        @. f.s[g.ie+1, :, :] = 2 * s_east - f.s[g.ie, :, :]

        interp = interpolate((gsrc.x, gsrc.y), fsrc.s_gradbot, (Gridded(Linear()), Gridded(Linear())))
        f.s_gradbot[:, :] .= interp(g.x, g.y)
        interp = interpolate((gsrc.x, gsrc.y), fsrc.s_gradtop, (Gridded(Linear()), Gridded(Linear())))
        f.s_gradtop[:, :] .= interp(g.x, g.y)

        # Set the pressure to that of the parent domain to correct only the pressure gradient.
        interp = interpolate((gsrc.x, gsrc.y, gsrc.z), fsrc.p, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
        f.p[:, :, :] .= interp(g.x, g.y, g.z)

    else
        g = m.grid[i]
        f = m.fields[i]

        # Set the walls to const normal velocity (zero for no-penetration BC).
        f.u[g.is  , :, :] .= 0
        f.u[g.ie+1, :, :] .= 0
        # z_tmp = reshape(g.z[:], (1, g.kcells))
        # zsize = g.zsize
        # f.u[g.is  , :, :] .= z_tmp ./ zsize
        # f.u[g.ie+1, :, :] .= z_tmp ./ zsize

        # For other velocity components we impose free slip (no-flux).
        f.v[g.is-1, :, :] .= f.v[g.is, :, :]
        f.v[g.ie+1, :, :] .= f.v[g.ie, :, :]
        f.w[g.is-1, :, :] .= f.w[g.is, :, :]
        f.w[g.ie+1, :, :] .= f.w[g.ie, :, :]
        f.s[g.is-1, :, :] .= f.s[g.is, :, :]
        f.s[g.ie+1, :, :] .= f.s[g.ie, :, :]

        # Set the BCs again, this needs to be done more clean if ever implemented properly...
        b = m.boundary[1]
        # Bottom BC.
        set_ghost_cells_bot_kernel!(
            f.u, f.u_bot, f.u_gradbot, g.dzh, g.ks, b.mom_bot_type)
        set_ghost_cells_bot_kernel!(
            f.v, f.v_bot, f.v_gradbot, g.dzh, g.ks, b.mom_bot_type)
        set_ghost_cells_bot_kernel!(
            f.s, f.s_bot, f.s_gradbot, g.dzh, g.ks, b.s_bot_type)

        # Top BC.
        set_ghost_cells_top_kernel!(
            f.u, f.u_top, f.u_gradtop, g.dzh, g.ke, b.mom_top_type)
        set_ghost_cells_top_kernel!(
            f.v, f.v_top, f.v_gradtop, g.dzh, g.ke, b.mom_top_type)
        set_ghost_cells_top_kernel!(
            f.s, f.s_top, f.s_gradtop, g.dzh, g.ke, b.s_top_type)
    end
    # CvH END TMP

    @timeit m.to "calc_dynamics_tend" calc_dynamics_tend!(m.fields[i], m.grid[i], m.to)
    @timeit m.to "calc_nudge_tend"    calc_nudge_tend!(m.fields[i], m.grid[i], m.multidomain[i])
    @timeit m.to "calc_pressure_tend" calc_pressure_tend!(m.fields[i], m.grid[i], m.timeloop[i], m.pressure[i], m.boundary[i], m.parallel, m.to)
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
        write(fid, "u_bot", f.u_bot[g.is:g.ie, g.js:g.je])
        write(fid, "u_top", f.u_top[g.is:g.ie, g.js:g.je])
        write(fid, "u_gradbot", f.u_gradbot[g.is:g.ie, g.js:g.je])
        write(fid, "u_gradtop", f.u_gradtop[g.is:g.ie, g.js:g.je])
        write(fid, "v_bot", f.v_bot[g.is:g.ie, g.js:g.je])
        write(fid, "v_top", f.v_top[g.is:g.ie, g.js:g.je])
        write(fid, "v_gradbot", f.v_gradbot[g.is:g.ie, g.js:g.je])
        write(fid, "v_gradtop", f.v_gradtop[g.is:g.ie, g.js:g.je])
        write(fid, "s_bot", f.s_bot[g.is:g.ie, g.js:g.je])
        write(fid, "s_top", f.s_top[g.is:g.ie, g.js:g.je])
        write(fid, "s_gradbot", f.s_gradbot[g.is:g.ie, g.js:g.je])
        write(fid, "s_gradtop", f.s_gradtop[g.is:g.ie, g.js:g.je])
        write(fid, "s_ref", f.s_ref[g.ks:g.ke])

        # Attach the dimensions. Note the c-indexing,
        # zero-based and swapped dim order.
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

        HDF5.h5ds_attach_scale(fid["u_bot"], fid["xh"], 1)
        HDF5.h5ds_attach_scale(fid["u_bot"], fid["y"], 0)

        HDF5.h5ds_attach_scale(fid["u_top"], fid["xh"], 1)
        HDF5.h5ds_attach_scale(fid["u_top"], fid["y"], 0)

        HDF5.h5ds_attach_scale(fid["u_gradbot"], fid["xh"], 1)
        HDF5.h5ds_attach_scale(fid["u_gradbot"], fid["y"], 0)

        HDF5.h5ds_attach_scale(fid["u_gradtop"], fid["xh"], 1)
        HDF5.h5ds_attach_scale(fid["u_gradtop"], fid["y"], 0)

        HDF5.h5ds_attach_scale(fid["v_bot"], fid["x"], 1)
        HDF5.h5ds_attach_scale(fid["v_bot"], fid["yh"], 0)

        HDF5.h5ds_attach_scale(fid["v_top"], fid["x"], 1)
        HDF5.h5ds_attach_scale(fid["v_top"], fid["yh"], 0)

        HDF5.h5ds_attach_scale(fid["v_gradbot"], fid["x"], 1)
        HDF5.h5ds_attach_scale(fid["v_gradbot"], fid["yh"], 0)

        HDF5.h5ds_attach_scale(fid["v_gradtop"], fid["x"], 1)
        HDF5.h5ds_attach_scale(fid["v_gradtop"], fid["yh"], 0)

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

        # Remove the ghost cells.
        a_nogc = a[g.is:g.ie, g.js:g.je, g.ks:g.ks+ktot-1]

        # Create the dataset for the entire field.
        aid = create_dataset(fid, name, datatype(a_nogc), dataspace((g.itot, g.jtot, ktot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)

        # Save the top slice before running transposes.
        if ktot == g.ktoth
            aid[is:ie, js:je, ktot] = a_nogc[:, :, ktot]
        end

        # Take a view without top.
        a_nogc_notop = @view a_nogc[:, :, 1:g.ktot]

        # Transpose the data to xz format to have far more contiguous data.
        # This prevents having to use chunking and keeps files viewable.
        a_nogc_t = reshape(a_nogc_notop, (g.itot, g.jmax, g.kblock))
        transpose_zx!(a_nogc_t, a_nogc_notop, g, p)

        # Save the data to disk.
        ks = p.id_x*g.kblock + 1; ke = (p.id_x+1)*g.kblock
        aid[:, js:je, ks:ke] = a_nogc_t[:, :, 1:g.kblock]

        close(aid)
    end

    items_to_save = [
        ("u_bot", f.u_bot),
        ("u_top", f.u_top),
        ("u_gradbot", f.u_gradbot),
        ("u_gradtop", f.u_gradtop),
        ("v_bot", f.v_bot),
        ("v_top", f.v_top),
        ("v_gradbot", f.v_gradbot),
        ("v_gradtop", f.v_gradtop),
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
        ("u_bot", "xh", "y"),
        ("u_top", "xh", "y"),
        ("u_gradbot", "xh", "y"),
        ("u_gradtop", "xh", "y"),
        ("v_bot", "x", "yh"),
        ("v_top", "x", "yh"),
        ("v_gradbot", "x", "yh"),
        ("v_gradtop", "x", "yh"),
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
        f.u_bot[g.is:g.ie, g.js:g.je] = read(fid, "u_bot")
        f.u_top[g.is:g.ie, g.js:g.je] = read(fid, "u_top")
        f.u_gradbot[g.is:g.ie, g.js:g.je] = read(fid, "u_gradbot")
        f.u_gradtop[g.is:g.ie, g.js:g.je] = read(fid, "u_gradtop")
        f.v_bot[g.is:g.ie, g.js:g.je] = read(fid, "v_bot")
        f.v_top[g.is:g.ie, g.js:g.je] = read(fid, "v_top")
        f.v_gradbot[g.is:g.ie, g.js:g.je] = read(fid, "v_gradbot")
        f.v_gradtop[g.is:g.ie, g.js:g.je] = read(fid, "v_gradtop")
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

    u_bot_id = open_dataset(fid, "u_bot", dapl, dxpl)
    f.u_bot[g.is:g.ie, g.js:g.je] = u_bot_id[is:ie, js:je]
    close(s_bot_id)

    u_top_id = open_dataset(fid, "u_top", dapl, dxpl)
    f.u_top[g.is:g.ie, g.js:g.je] = u_top_id[is:ie, js:je]
    close(u_top_id)

    u_gradbot_id = open_dataset(fid, "u_gradbot", dapl, dxpl)
    f.u_gradbot[g.is:g.ie, g.js:g.je] = u_gradbot_id[is:ie, js:je]
    close(u_gradbot_id)

    u_gradtop_id = open_dataset(fid, "u_gradtop", dapl, dxpl)
    f.u_gradtop[g.is:g.ie, g.js:g.je] = u_gradtop_id[is:ie, js:je]
    close(u_gradtop_id)

    v_bot_id = open_dataset(fid, "v_bot", dapl, dxpl)
    f.v_bot[g.is:g.ie, g.js:g.je] = v_bot_id[is:ie, js:je]
    close(v_bot_id)

    v_top_id = open_dataset(fid, "v_top", dapl, dxpl)
    f.v_top[g.is:g.ie, g.js:g.je] = v_top_id[is:ie, js:je]
    close(v_top_id)

    v_gradbot_id = open_dataset(fid, "v_gradbot", dapl, dxpl)
    f.v_gradbot[g.is:g.ie, g.js:g.je] = v_gradbot_id[is:ie, js:je]
    close(v_gradbot_id)

    v_gradtop_id = open_dataset(fid, "v_gradtop", dapl, dxpl)
    f.v_gradtop[g.is:g.ie, g.js:g.je] = v_gradtop_id[is:ie, js:je]
    close(v_gradtop_id)

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

    # CvH TMP disable zoom-in looping
    # @timeit m.to "domains" begin
    #     for i in 1:m.n_domains
    #         # Compute the nudging fields, which remain constant throughout the step.
    #         if i > 1
    #             @timeit m.to "calc_nudge_fields" calc_nudge_fields!(m.multidomain[i], m.fields[i], m.fields[i-1], m.grid[i], m.grid[i-1])
    #         end

    #         while (m.timeloop[i].itime < itime_next)
    #             @timeit m.to "integrate_time" integrate_time!(m.fields[i], m.grid[i], m.timeloop[i])
    #             @timeit m.to "step_time" step_time!(m.timeloop[i])

    #             if is_save_time(m.timeloop[i])
    #                 @timeit m.to "save_domain" save_domain(m, i, m.parallel)
    #             end

    #             @timeit m.to "calc_rhs" calc_rhs!(m, i)
    #         end
    #     end
    # end

    # if is_check_time(m.timeloop[1])
    #     @timeit m.to "check_model" check_model(m)
    # end
    # CvH TMP disable zoom-in looping

    # CvH TMP
    # Step all domains with the same time step.
    @timeit m.to "domains" begin
        while (m.timeloop[1].itime < itime_next)
            for i in 1:m.n_domains
                @timeit m.to "integrate_time" integrate_time!(m.fields[i], m.grid[i], m.timeloop[i])
            end

            # CvH TMP TWO WAY NEST HERE

            if m.n_domains > 1
                f1 = m.fields[1]; g1 = m.grid[1]
                f2 = m.fields[2]; g2 = m.grid[2]

                # HARDCODE THE NESTING, FIX THIS ONCE IT WORKS
                isc = g1.is + g1.itot÷4; iec = g1.is + 3*(g1.itot÷4) - 1
                jsc = g1.js; jec = g1.je
                ksc = g1.ks; kec = g1.ke

                ni = round(Int, g1.dx / g2.dx)
                nj = round(Int, g1.dy / g2.dy)
                nk = round(Int, g1.dz[g1.ks] / g2.dz[g2.ks])

                for k in ksc:kec
                    for j in jsc:jec
                        for i in isc:iec
                            iis = 2*(i-isc) + g2.is; jjs = 2*(j-jsc) + g2.js; kks = 2*(k-ksc) + g2.ks
                            iie = iis+ni-1; jje = jjs+nj-1; kke = kks+nk-1;
                            f1.s[i, j, k] = sum(f2.s[iis:iie, jjs:jje, kks:kke]) / (ni*nj*nk);
                        end
                    end
                end

                # for k in ksc:kec
                #     for j in jsc:jec
                #         for i in isc:iec
                #             ii = 2*(i-isc); jj = 2*(j-jsc); kk = 2*(k-ksc)
                #             f1.u[i, j, k] = 0;
                #         end
                #     end
                # end

                # for k in ksc:kec
                #     for j in jsc:jec
                #         for i in isc:iec
                #             ii = 2*(i-isc); jj = 2*(j-jsc); kk = 2*(k-ksc)
                #             f1.v[i, j, k] = 0;
                #         end
                #     end
                # end

                # for k in ksc:kec
                #     for j in jsc:jec
                #         for i in isc:iec
                #             ii = 2*(i-isc); jj = 2*(j-jsc); kk = 2*(k-ksc)
                #             f1.w[i, j, k] = 0;
                #         end
                #     end
                # end

                # interp = interpolate((g2.xh, g2.y, g2.z), f2.u, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
                # f1.u[isc:iec+1, jsc:jec, ksc:kec] .= interp(g1.xh[isc:iec+1], g1.y[jsc:jec], g1.z[ksc:kec])

                # interp = interpolate((g2.x, g2.yh, g2.z), f2.v, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
                # f1.v[isc:iec, jsc:jec, ksc:kec] .= interp(g1.x[isc:iec], g1.yh[jsc:jec], g1.z[ksc:kec])

                # interp = interpolate((g2.x, g2.y, g2.zh), f2.w, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
                # f1.w[isc:iec, jsc:jec, ksc:kec] .= interp(g1.x[isc:iec], g1.y[jsc:jec], g1.zh[ksc:kec])

                # interp = interpolate((g2.x, g2.y, g2.z), f2.s, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
                # f1.s[isc:iec, jsc:jec, ksc:kec] .= interp(g1.x[isc:iec], g1.y[jsc:jec], g1.z[ksc:kec])
            end
            # CvH END TMP

            for i in 1:m.n_domains
                @timeit m.to "step_time" step_time!(m.timeloop[i])
            end

            for i in 1:m.n_domains
                if is_save_time(m.timeloop[i])
                    @timeit m.to "save_domain" save_domain(m, i, m.parallel)
                end
            end

            for i in 1:m.n_domains
                @timeit m.to "calc_rhs" calc_rhs!(m, i)
            end
        end
    end

    if is_check_time(m.timeloop[1])
        @timeit m.to "check_model" check_model(m)
    end
    # CvH TMP END

    # Enable the timer after the first round to avoid measuring compilation.
    enable_timer!(m.to)

    return is_in_progress(m.timeloop[1])
end


function check_model(m::Model)
    # Calculate the time since the last check.
    old_time = m.last_measured_time[]
    m.last_measured_time[] = time_ns()

    # First, print the model time and the wall clock since last sample.
    status_string = @sprintf("(%11.3f) Time = %8.3f",
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
        println("")
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
