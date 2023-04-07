struct Fields{TF <: Union{Float32, Float64}}
    u::Array{TF, 3}
    v::Array{TF, 3}
    w::Array{TF, 3}
    s::Array{TF, 3}
    scalars::Dict{String, Array{TF, 3}}

    u_tend::Array{TF, 3}
    v_tend::Array{TF, 3}
    w_tend::Array{TF, 3}
    s_tend::Array{TF, 3}
    scalars_tend::Dict{String, Array{TF, 3}}

    u_bot::Array{TF, 2}
    v_bot::Array{TF, 2}
    w_bot::Array{TF, 2}
    s_bot::Array{TF, 2}
    scalars_bot::Dict{String, Array{TF, 2}}

    u_gradbot::Array{TF, 2}
    v_gradbot::Array{TF, 2}
    w_gradbot::Array{TF, 2}
    s_gradbot::Array{TF, 2}
    scalars_gradbot::Dict{String, Array{TF, 2}}

    u_top::Array{TF, 2}
    v_top::Array{TF, 2}
    w_top::Array{TF, 2}
    s_top::Array{TF, 2}
    scalars_top::Dict{String, Array{TF, 2}}

    u_gradtop::Array{TF, 2}
    v_gradtop::Array{TF, 2}
    w_gradtop::Array{TF, 2}
    s_gradtop::Array{TF, 2}
    scalars_gradtop::Dict{String, Array{TF, 2}}

    s_ref::Array{TF, 1}

    p::Array{TF, 3}

    visc::TF
    alpha::TF
end


function Fields(g::Grid, settings::Dict, TF)
    d = settings["fields"]

    visc = d["visc"]
    alpha = d["alpha"]

    # Read in passive scalars.
    scalars = Dict{String, Array{TF, 3}}()
    scalars_tend = Dict{String, Array{TF, 3}}()
    scalars_bot = Dict{String, Array{TF, 2}}()
    scalars_gradbot = Dict{String, Array{TF, 2}}()
    scalars_top = Dict{String, Array{TF, 2}}()
    scalars_gradtop = Dict{String, Array{TF, 2}}()

    if haskey(d, "scalars")
        scalar_names = d["scalars"]

        for scalar_name in scalar_names
            scalars[scalar_name] = zeros(g.icells, g.jcells, g.kcells)
            scalars_tend[scalar_name] = zeros(g.icells, g.jcells, g.kcells)
            scalars_bot[scalar_name] = zeros(g.icells, g.jcells)
            scalars_top[scalar_name] = zeros(g.icells, g.jcells)
            scalars_gradbot[scalar_name] = zeros(g.icells, g.jcells)
            scalars_gradtop[scalar_name] = zeros(g.icells, g.jcells)
        end
    end
    # End of passive scalar processing.


    Fields{TF}(
        zeros(g.icells, g.jcells, g.kcells), # u
        zeros(g.icells, g.jcells, g.kcells), # v
        zeros(g.icells, g.jcells, g.kcells), # w
        zeros(g.icells, g.jcells, g.kcells), # s
        scalars,

        zeros(g.icells, g.jcells, g.kcells), # u_tend
        zeros(g.icells, g.jcells, g.kcells), # v_tend
        zeros(g.icells, g.jcells, g.kcells), # w_tend
        zeros(g.icells, g.jcells, g.kcells), # s_tend
        scalars_tend,

        zeros(g.icells, g.jcells), # u_bot
        zeros(g.icells, g.jcells), # v_bot
        zeros(g.icells, g.jcells), # w_bot
        zeros(g.icells, g.jcells), # s_bot
        scalars_bot,

        zeros(g.icells, g.jcells), # u_gradbot
        zeros(g.icells, g.jcells), # v_gradbot
        zeros(g.icells, g.jcells), # w_gradbot
        zeros(g.icells, g.jcells), # s_gradbot
        scalars_gradbot,

        zeros(g.icells, g.jcells), # u_top
        zeros(g.icells, g.jcells), # v_top
        zeros(g.icells, g.jcells), # w_top
        zeros(g.icells, g.jcells), # s_top
        scalars_top,

        zeros(g.icells, g.jcells), # u_gradtop
        zeros(g.icells, g.jcells), # v_gradtop
        zeros(g.icells, g.jcells), # w_gradtop
        zeros(g.icells, g.jcells), # s_gradtop
        scalars_gradtop,

        zeros(g.kcells), # s_ref

        zeros(g.icells, g.jcells, g.kcells), # p

        visc, alpha)
end
