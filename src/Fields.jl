struct Fields
    u::Array{Float64, 3}
    v::Array{Float64, 3}
    w::Array{Float64, 3}
    s::Array{Float64, 3}

    u_tend::Array{Float64, 3}
    v_tend::Array{Float64, 3}
    w_tend::Array{Float64, 3}
    s_tend::Array{Float64, 3}

    p::Array{Float64, 3}

    visc::Float64
end

function Fields(g::Grid, d::Dict)
    visc = d["visc"]

    Fields(
        zeros(g.icells, g.jcells, g.kcells),
        zeros(g.icells, g.jcells, g.kcells),
        zeros(g.icells, g.jcells, g.kcells),
        zeros(g.icells, g.jcells, g.kcells),

        zeros(g.icells, g.jcells, g.kcells),
        zeros(g.icells, g.jcells, g.kcells),
        zeros(g.icells, g.jcells, g.kcells),
        zeros(g.icells, g.jcells, g.kcells),

        zeros(g.icells, g.jcells, g.kcells),

        visc)
end
