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

function set_ghost_cells_bot_kernel!(
    a, a_bot, a_gradbot, dz, ks, bot_type::Boundary_type)

    if bot_type == Dirichlet::Boundary_type
        @inbounds @. a[:, :, ks-1] = 2a_bot[:, :] - a[:, :, ks]
    elseif bot_type == Neumann::Boundary_type
        @inbounds @. a[:, :, ks-1] = a_gradbot[:, :] + a[:, :, ks]
    end
end

function set_ghost_cells_top_kernel!(
    a, a_top, a_gradtop, dz, ke, bot_type::Boundary_type)

    if bot_type == Dirichlet::Boundary_type
        @inbounds @. a[:, :, ke+1] = 2a_top[:, :] - a[:, :, ke]
    elseif bot_type == Neumann::Boundary_type
        @inbounds @. a[:, :, ke+1] = a_gradtop[:, :]*dz + a[:, :, ke]
    end
end

function set_boundary!(fields::Fields, grid::Grid, boundary::Boundary)
    # Bottom BC.
    set_ghost_cells_bot_kernel!(
        fields.u, fields.u_bot, fields.u_gradbot, grid.dz, grid.ks, boundary.mom_bot_type)
    set_ghost_cells_bot_kernel!(
        fields.v, fields.v_bot, fields.v_gradbot, grid.dz, grid.ks, boundary.mom_bot_type)
    set_ghost_cells_bot_kernel!(
        fields.s, fields.s_bot, fields.s_gradbot, grid.dz, grid.ks, boundary.s_bot_type)

    # Top BC.
    set_ghost_cells_top_kernel!(
        fields.u, fields.u_top, fields.u_gradtop, grid.dz, grid.ke, boundary.mom_top_type)
    set_ghost_cells_top_kernel!(
        fields.v, fields.v_top, fields.v_gradtop, grid.dz, grid.ke, boundary.mom_top_type)
    set_ghost_cells_top_kernel!(
        fields.s, fields.s_top, fields.s_gradtop, grid.dz, grid.ke, boundary.s_top_type)
end
