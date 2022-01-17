@enum Boundary_type Dirichlet Neumann

struct Boundary
    mom_bot_type::Boundary_type
    mom_top_type::Boundary_type
    s_bot_type::Boundary_type
    s_top_type::Boundary_type
end


function string_to_type(s::String)
    if s == "Dirichlet"
        return Dirichlet::Boundary_type
    elseif s == "Neumann"
        return Neumann::Boundary_type
    end
end


function Boundary(d::Dict)
    mom_bot_type = string_to_type(d["mom_bot_type"])
    mom_top_type = string_to_type(d["mom_top_type"])
    s_bot_type = string_to_type(d["s_bot_type"])
    s_top_type = string_to_type(d["s_top_type"])

    Boundary(
        mom_bot_type,
        mom_top_type,
        s_bot_type,
        s_top_type)
end


function boundary_cyclic_kernel!(
    a, is, ie, js, je, igc, jgc, p::ParallelSerial)
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
    a, is, ie, js, je, igc, jgc, p::ParallelDistributed)

    # Transfer the east-west data.
    send_east = a[ie-igc+1:ie, :, :]
    send_west = a[is:is+igc-1, :, :]
    recv_west = similar(send_east)
    recv_east = similar(send_west)

    reqs = Vector{MPI.Request}(undef, 4)
    reqs[1] = MPI.Irecv!(recv_west, p.id_west, 1, p.commxy);
    reqs[2] = MPI.Irecv!(recv_east, p.id_east, 2, p.commxy);
    reqs[3] = MPI.Isend( send_east, p.id_east, 1, p.commxy);
    reqs[4] = MPI.Isend( send_west, p.id_west, 2, p.commxy);
    MPI.Waitall!(reqs)

    @tturbo a[   1:igc, :, :] .= recv_west[:, :, :]
    @tturbo a[ie+1:end, :, :] .= recv_east[:, :, :]

    # Transfer the north-south data.
    send_north = a[:, je-jgc+1:je, :]
    send_south = a[:, js:js+jgc-1, :]
    recv_south = similar(send_north)
    recv_north = similar(send_south)

    reqs[1] = MPI.Irecv!(recv_north, p.id_north, 1, p.commxy);
    reqs[2] = MPI.Irecv!(recv_south, p.id_south, 2, p.commxy);
    reqs[3] = MPI.Isend( send_south, p.id_south, 1, p.commxy);
    reqs[4] = MPI.Isend( send_north, p.id_north, 2, p.commxy);
    MPI.Waitall!(reqs)

    @tturbo a[:,    1:jgc, :] .= recv_south[:, :, :]
    @tturbo a[:, je+1:end, :] .= recv_north[:, :, :]
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
    a, a_top, a_gradtop, dzh, ke, bot_type::Boundary_type)

    if bot_type == Dirichlet::Boundary_type
        dzhi = 1 / dzh[ke+1]
        @tturbo for j in 1:size(a, 2)
            for i in 1:size(a, 1)
                a[i, j, ke+1] = 2a_top[i, j] - a[i, j, ke]
                a_gradtop[i, j] = (a[i, j, ke+1] - a[i, j, ke]) * dzhi
            end
        end

    elseif bot_type == Neumann::Boundary_type
        @tturbo for j in 1:size(a, 2)
            for i in 1:size(a, 1)
                a[i, j, ke+1] = a_gradtop[i, j]*dzh[ke+1] + a[i, j, ke]
                a_top[i, j] = 1//2 * (a[i, j, ke] + a[i, j, ke+1])
            end
        end
    end
end


function set_boundary!(fields::Fields, grid::Grid, boundary::Boundary, p::Parallel)
    # Cyclic BC.
    boundary_cyclic_kernel!(
        fields.u, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.v, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.w, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.s, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)

    # Cyclic BCs of boundary fields.
    boundary_cyclic_kernel!(
        fields.u_bot, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.u_gradbot, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.u_top, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.u_gradtop, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)

    boundary_cyclic_kernel!(
        fields.v_bot, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.v_gradbot, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.v_top, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.v_gradtop, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)

    boundary_cyclic_kernel!(
        fields.s_bot, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.s_gradbot, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.s_top, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    boundary_cyclic_kernel!(
        fields.s_gradtop, grid.is, grid.ie, grid.js, grid.je, grid.igc, grid.jgc, p)
    
    # Bottom BC.
    set_ghost_cells_bot_kernel!(
        fields.u, fields.u_bot, fields.u_gradbot, grid.dzh, grid.ks, boundary.mom_bot_type)
    set_ghost_cells_bot_kernel!(
        fields.v, fields.v_bot, fields.v_gradbot, grid.dzh, grid.ks, boundary.mom_bot_type)
    set_ghost_cells_bot_kernel!(
        fields.s, fields.s_bot, fields.s_gradbot, grid.dzh, grid.ks, boundary.s_bot_type)

    # Top BC.
    set_ghost_cells_top_kernel!(
        fields.u, fields.u_top, fields.u_gradtop, grid.dzh, grid.ke, boundary.mom_top_type)
    set_ghost_cells_top_kernel!(
        fields.v, fields.v_top, fields.v_gradtop, grid.dzh, grid.ke, boundary.mom_top_type)
    set_ghost_cells_top_kernel!(
        fields.s, fields.s_top, fields.s_gradtop, grid.dzh, grid.ke, boundary.s_top_type)
end
