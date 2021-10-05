using Interpolations

struct MultiDomain{TF <: Union{Float32, Float64}}
    enable_nudge::Bool
    nudge_time::TF

    u_nudge::Array{TF, 3}
    v_nudge::Array{TF, 3}
    w_nudge::Array{TF, 3}
    s_nudge::Array{TF, 3}
end


# Constructor for multidomain struct.
function MultiDomain(g::Grid, settings, TF)
    enable_nudge = settings["enable_nudge"]
    if enable_nudge
        nudge_time = settings["nudge_time"]
        return MultiDomain{TF}(
            enable_nudge,
            nudge_time,
            zeros(g.icells, g.jcells, g.kcells),
            zeros(g.icells, g.jcells, g.kcells),
            zeros(g.icells, g.jcells, g.kcells),
            zeros(g.icells, g.jcells, g.kcells))
    else
        nudge_time = 1e12
        return MultiDomain{TF}(
            enable_nudge,
            nudge_time,
            zeros(0, 0, 0),
            zeros(0, 0, 0),
            zeros(0, 0, 0),
            zeros(0, 0, 0))
    end
end


# Transfer the state
function transfer_state!(f_d::Fields, f_s::Fields, g_d::Grid, g_s::Grid)
    interp = LinearInterpolation((g_s.xh, g_s.y, g_s.z), f_s.u)
    f_d.u[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] = interp(g_d.xh[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = LinearInterpolation((g_s.x, g_s.yh, g_s.z), f_s.v)
    f_d.v[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] = interp(g_d.x[g_d.is:g_d.ie], g_d.yh[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = LinearInterpolation((g_s.x, g_s.y, g_s.zh), f_s.w)
    f_d.w[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.keh] = interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.zh[g_d.ks:g_d.keh])

    interp = LinearInterpolation((g_s.x, g_s.y, g_s.z), f_s.s)
    f_d.s[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] = interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = LinearInterpolation((g_s.x, g_s.y, g_s.z), f_s.p)
    f_d.p[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] = interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = LinearInterpolation((g_s.x, g_s.y), f_s.s_bot)
    f_d.s_bot[g_d.is:g_d.ie, g_d.js:g_d.je] = interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je])

    interp = LinearInterpolation((g_s.x, g_s.y), f_s.s_gradbot)
    f_d.s_gradbot[g_d.is:g_d.ie, g_d.js:g_d.je] = interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je])
end


function calc_nudge_fields!(md::MultiDomain, f_d::Fields, f_s::Fields, g_d::Grid, g_s::Grid)
    if !md.enable_nudge
        return
    end

    interp = LinearInterpolation((g_s.xh, g_s.y, g_s.z), f_s.u)
    md.u_nudge[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] = interp(g_d.xh[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = LinearInterpolation((g_s.x, g_s.yh, g_s.z), f_s.v)
    md.v_nudge[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] = interp(g_d.x[g_d.is:g_d.ie], g_d.yh[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = LinearInterpolation((g_s.x, g_s.y, g_s.zh), f_s.w)
    md.w_nudge[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.keh] = interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.zh[g_d.ks:g_d.keh])

    interp = LinearInterpolation((g_s.x, g_s.y, g_s.z), f_s.s)
    md.s_nudge[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] = interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])
end


function nudge_kernel!(a_tend, a, a_nudge, c_nudge, is, ie, js, je, ks, ke)
    @fast3d begin
        @fd (a_tend, a, a_nudge) a_tend -= c_nudge * (a - a_nudge)
    end
end


function calc_nudge_tend!(f::Fields, g::Grid, md::MultiDomain)
    if !md.enable_nudge
        return
    end

    c_nudge = 1 / md.nudge_time

    nudge_kernel!(f.u_tend, f.u, md.u_nudge, c_nudge, g.is, g.ie, g.js, g.je, g.ks, g.ke )
    nudge_kernel!(f.v_tend, f.v, md.v_nudge, c_nudge, g.is, g.ie, g.js, g.je, g.ks, g.ke )
    nudge_kernel!(f.w_tend, f.w, md.w_nudge, c_nudge, g.is, g.ie, g.js, g.je, g.ks, g.keh)
    nudge_kernel!(f.s_tend, f.s, md.s_nudge, c_nudge, g.is, g.ie, g.js, g.je, g.ks, g.ke )
end
