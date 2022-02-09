struct Fields{TF <: Union{Float32, Float64}}
    u::Array{TF, 3}
    v::Array{TF, 3}
    w::Array{TF, 3}
    s::Array{TF, 3}

    u_tend::Array{TF, 3}
    v_tend::Array{TF, 3}
    w_tend::Array{TF, 3}
    s_tend::Array{TF, 3}

    u_bot::Array{TF, 2}
    v_bot::Array{TF, 2}
    w_bot::Array{TF, 2}
    s_bot::Array{TF, 2}

    u_gradbot::Array{TF, 2}
    v_gradbot::Array{TF, 2}
    w_gradbot::Array{TF, 2}
    s_gradbot::Array{TF, 2}

    u_top::Array{TF, 2}
    v_top::Array{TF, 2}
    w_top::Array{TF, 2}
    s_top::Array{TF, 2}

    u_gradtop::Array{TF, 2}
    v_gradtop::Array{TF, 2}
    w_gradtop::Array{TF, 2}
    s_gradtop::Array{TF, 2}

    s_ref::Array{TF, 1}

    p::Array{TF, 3}

    visc::TF
    alpha::TF
end


function Fields(g::Grid, d::Dict, TF)
    visc = d["visc"]
    alpha = d["alpha"]

    Fields{TF}(
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

        zeros(g.kcells), # s_ref

        zeros(g.icells, g.jcells, g.kcells), # p

        visc, alpha)
end
