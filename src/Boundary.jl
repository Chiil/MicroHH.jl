@enum Boundary_type Dirichlet Neumann

struct BoundaryBuffers{TF <: Union{Float32, Float64}}
    send_west::Array{TF, 3}
    recv_west::Array{TF, 3}
    send_east::Array{TF, 3}
    recv_east::Array{TF, 3}

    send_south::Array{TF, 3}
    recv_south::Array{TF, 3}
    send_north::Array{TF, 3}
    recv_north::Array{TF, 3}
end


struct Boundary
    mom_bot_type::Boundary_type
    mom_top_type::Boundary_type
    s_bot_type::Boundary_type
    s_top_type::Boundary_type

    buffers::BoundaryBuffers
end


function string_to_type(s::String)
    if s == "Dirichlet"
        return Dirichlet::Boundary_type
    elseif s == "Neumann"
        return Neumann::Boundary_type
    end
end


function Boundary(g::Grid, pp::Parallel, settings::Dict, TF)
    d = settings["boundary"]

    mom_bot_type = string_to_type(d["mom_bot_type"])
    mom_top_type = string_to_type(d["mom_top_type"])
    s_bot_type = string_to_type(d["s_bot_type"])
    s_top_type = string_to_type(d["s_top_type"])

    if pp isa ParallelSerial
        buffers = BoundaryBuffers(
            zeros(TF, (0, 0, 0)),
            zeros(TF, (0, 0, 0)),
            zeros(TF, (0, 0, 0)),
            zeros(TF, (0, 0, 0)),

            zeros(TF, (0, 0, 0)),
            zeros(TF, (0, 0, 0)),
            zeros(TF, (0, 0, 0)),
            zeros(TF, (0, 0, 0)) )
    else pp isa ParallelDistributed
        buffers = BoundaryBuffers(
            zeros(TF, (g.igc, g.jcells, g.kcells)),
            zeros(TF, (g.igc, g.jcells, g.kcells)),
            zeros(TF, (g.igc, g.jcells, g.kcells)),
            zeros(TF, (g.igc, g.jcells, g.kcells)),

            zeros(TF, (g.icells, g.jgc, g.kcells)),
            zeros(TF, (g.icells, g.jgc, g.kcells)),
            zeros(TF, (g.icells, g.jgc, g.kcells)),
            zeros(TF, (g.icells, g.jgc, g.kcells)) )
    end

    Boundary(
        mom_bot_type, mom_top_type,
        s_bot_type, s_top_type,
        buffers)
end


function boundary_cyclic_kernel!(
    a, is, ie, js, je, igc, jgc, bufs::BoundaryBuffers, p::ParallelSerial)

    # East-west BCs
    @tturbo for k in 1:size(a, 3)
        for j in 1:size(a, 2)
            for i in 1:igc
                a[i, j, k] = a[ie-i+1, j, k]
                a[ie+i, j, k] = a[is+i-1, j, k]
            end
        end
    end

    # North-south BCs
    @tturbo for k in 1:size(a, 3)
        for j in 1:jgc
            for i in 1:size(a, 1)
                a[i, j, k] = a[i, je-j+1, k]
                a[i, je+j, k] = a[i, js+j-1, k]
            end
        end
    end
end


function boundary_cyclic_kernel!(
    a, is, ie, js, je, igc, jgc, buf::BoundaryBuffers, p::ParallelDistributed)

    if length(size(a)) == 2 || (length(size(a)) == 3 && size(a, 3) == 1)
        sw = @view buf.send_west[:, :, 1:1]
        rw = @view buf.recv_west[:, :, 1:1]
        se = @view buf.send_east[:, :, 1:1]
        re = @view buf.recv_east[:, :, 1:1]
        ss = @view buf.send_south[:, :, 1:1]
        rs = @view buf.recv_south[:, :, 1:1]
        sn = @view buf.send_north[:, :, 1:1]
        rn = @view buf.recv_north[:, :, 1:1]
    else
        sw = buf.send_west; rw = buf.recv_west
        se = buf.send_east; re = buf.recv_east
        ss = buf.send_south; rs = buf.recv_south
        sn = buf.send_north; rn = buf.recv_north
    end

    # Transfer the east-west data.
    a_east = @view a[ie-igc+1:ie, :, :]
    a_west = @view a[is:is+igc-1, :, :]
    @tturbo se .= a_east
    @tturbo sw .= a_west

    reqs = Vector{MPI.Request}(undef, 4)
    reqs[1] = MPI.Irecv!(rw, p.id_west, 1, p.commxy);
    reqs[2] = MPI.Irecv!(re, p.id_east, 2, p.commxy);
    reqs[3] = MPI.Isend( se, p.id_east, 1, p.commxy);
    reqs[4] = MPI.Isend( sw, p.id_west, 2, p.commxy);
    MPI.Waitall!(reqs)

    @tturbo a[   1:igc, :, :] .= rw
    @tturbo a[ie+1:end, :, :] .= re

    # Transfer the north-south data.
    a_north = @view a[:, je-jgc+1:je, :]
    a_south = @view a[:, js:js+jgc-1, :]
    @tturbo sn .= a_north
    @tturbo ss .= a_south

    reqs[1] = MPI.Irecv!(rn, p.id_north, 1, p.commxy);
    reqs[2] = MPI.Irecv!(rs, p.id_south, 2, p.commxy);
    reqs[3] = MPI.Isend( ss, p.id_south, 1, p.commxy);
    reqs[4] = MPI.Isend( sn, p.id_north, 2, p.commxy);
    MPI.Waitall!(reqs)

    @tturbo a[:,    1:jgc, :] .= rs
    @tturbo a[:, je+1:end, :] .= rn
end


function set_ghost_cells_bot_kernel!(
    a, a_bot, a_gradbot, dzh, ks, bot_type::Boundary_type)

    if bot_type == Dirichlet::Boundary_type
        dzhi = 1 / dzh[ks]
        @tturbo for j in 1:size(a, 2)
            for i in 1:size(a, 1)
                a[i, j, ks-1] = 2a_bot[i, j] - a[i, j, ks]
                a_gradbot[i, j] = (a[i, j, ks] - a[i, j, ks-1]) * dzhi
            end
        end

    elseif bot_type == Neumann::Boundary_type
        @tturbo for j in 1:size(a, 2)
            for i in 1:size(a, 1)
                a[i, j, ks-1] = -a_gradbot[i, j]*dzh[ks] + a[i, j, ks]
                a_bot[i, j] = 1//2 * (a[i, j, ks-1] + a[i, j, ks])
            end
        end
    end
end


function set_ghost_cells_top_kernel!(
    a, a_top, a_gradtop, dzh, ke, top_type::Boundary_type)

    if top_type == Dirichlet::Boundary_type
        dzhi = 1 / dzh[ke+1]
        @tturbo for j in 1:size(a, 2)
            for i in 1:size(a, 1)
                a[i, j, ke+1] = 2a_top[i, j] - a[i, j, ke]
                a_gradtop[i, j] = (a[i, j, ke+1] - a[i, j, ke]) * dzhi
            end
        end

    elseif top_type == Neumann::Boundary_type
        @tturbo for j in 1:size(a, 2)
            for i in 1:size(a, 1)
                a[i, j, ke+1] = a_gradtop[i, j]*dzh[ke+1] + a[i, j, ke]
                a_top[i, j] = 1//2 * (a[i, j, ke] + a[i, j, ke+1])
            end
        end
    end
end


function set_boundary!(f::Fields, g::Grid, b::Boundary, p::Parallel)
    # Cyclic BC.
    boundary_cyclic_kernel!(
        f.u, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.v, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.w, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.s, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)

    for scalar_name in keys(f.scalars)
        boundary_cyclic_kernel!(
            f.scalars[scalar_name], g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    end

    # Cyclic BCs of boundary f.
    boundary_cyclic_kernel!(
        f.u_bot, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.u_gradbot, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.u_top, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.u_gradtop, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)

    boundary_cyclic_kernel!(
        f.v_bot, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.v_gradbot, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.v_top, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.v_gradtop, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)

    boundary_cyclic_kernel!(
        f.s_bot, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.s_gradbot, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.s_top, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    boundary_cyclic_kernel!(
        f.s_gradtop, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)

    for scalar_name in keys(f.scalars)
        boundary_cyclic_kernel!(
            f.scalars_bot[scalar_name], g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
        boundary_cyclic_kernel!(
            f.scalars_gradbot[scalar_name], g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
        boundary_cyclic_kernel!(
            f.scalars_top[scalar_name], g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
        boundary_cyclic_kernel!(
            f.scalars_gradtop[scalar_name], g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, p)
    end

    # Bottom BC.
    set_ghost_cells_bot_kernel!(
        f.u, f.u_bot, f.u_gradbot, g.dzh, g.ks, b.mom_bot_type)
    set_ghost_cells_bot_kernel!(
        f.v, f.v_bot, f.v_gradbot, g.dzh, g.ks, b.mom_bot_type)
    set_ghost_cells_bot_kernel!(
        f.s, f.s_bot, f.s_gradbot, g.dzh, g.ks, b.s_bot_type)

    for scalar_name in keys(f.scalars)
        set_ghost_cells_bot_kernel!(
            f.scalars[scalar_name], f.scalars_bot[scalar_name], f.scalars_gradbot[scalar_name], g.dzh, g.ks, b.s_bot_type)
    end

    # Top BC.
    set_ghost_cells_top_kernel!(
        f.u, f.u_top, f.u_gradtop, g.dzh, g.ke, b.mom_top_type)
    set_ghost_cells_top_kernel!(
        f.v, f.v_top, f.v_gradtop, g.dzh, g.ke, b.mom_top_type)
    set_ghost_cells_top_kernel!(
        f.s, f.s_top, f.s_gradtop, g.dzh, g.ke, b.s_top_type)

    for scalar_name in keys(f.scalars)
        set_ghost_cells_top_kernel!(
            f.scalars[scalar_name], f.scalars_top[scalar_name], f.scalars_gradtop[scalar_name], g.dzh, g.ke, b.s_top_type)
    end
end
