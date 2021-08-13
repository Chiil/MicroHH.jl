struct Fields
    u::Array{Float64, 3}
    v::Array{Float64, 3}
    w::Array{Float64, 3}
    s::Array{Float64, 3}

    u_tend::Array{Float64, 3}
    v_tend::Array{Float64, 3}
    w_tend::Array{Float64, 3}
    s_tend::Array{Float64, 3}

    u_bot::Array{Float64, 2}
    v_bot::Array{Float64, 2}
    w_bot::Array{Float64, 2}
    s_bot::Array{Float64, 2}

    u_gradbot::Array{Float64, 2}
    v_gradbot::Array{Float64, 2}
    w_gradbot::Array{Float64, 2}
    s_gradbot::Array{Float64, 2}

    u_top::Array{Float64, 2}
    v_top::Array{Float64, 2}
    w_top::Array{Float64, 2}
    s_top::Array{Float64, 2}

    u_gradtop::Array{Float64, 2}
    v_gradtop::Array{Float64, 2}
    w_gradtop::Array{Float64, 2}
    s_gradtop::Array{Float64, 2}

    p::Array{Float64, 3}

    visc::Float64
end

function Fields(g::Grid, d::Dict)
    visc = d["visc"]

    Fields(
        zeros(g.icells, g.jcells, g.kcells), # u
        zeros(g.icells, g.jcells, g.kcells), # v
        zeros(g.icells, g.jcells, g.kcells), # w
        zeros(g.icells, g.jcells, g.kcells), # s

        zeros(g.icells, g.jcells, g.kcells), # u_tend
        zeros(g.icells, g.jcells, g.kcells), # v_tend
        zeros(g.icells, g.jcells, g.kcells), # w_tend
        zeros(g.icells, g.jcells, g.kcells), # s_tend

        zeros(g.icells, g.jcells), # u_bot
        zeros(g.icells, g.jcells), # v_bot
        zeros(g.icells, g.jcells), # w_bot
        zeros(g.icells, g.jcells), # s_bot

        zeros(g.icells, g.jcells), # u_gradbot
        zeros(g.icells, g.jcells), # v_gradbot
        zeros(g.icells, g.jcells), # w_gradbot
        zeros(g.icells, g.jcells), # s_gradbot

        zeros(g.icells, g.jcells), # u_top
        zeros(g.icells, g.jcells), # v_top
        zeros(g.icells, g.jcells), # w_top
        zeros(g.icells, g.jcells), # s_top

        zeros(g.icells, g.jcells), # u_gradtop
        zeros(g.icells, g.jcells), # v_gradtop
        zeros(g.icells, g.jcells), # w_gradtop
        zeros(g.icells, g.jcells), # s_gradtop

        zeros(g.icells, g.jcells, g.kcells), # p

        visc)
end
