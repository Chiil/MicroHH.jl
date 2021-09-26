using Interpolations

struct MultiDomain{TF <: Union{Float32, Float64}}
    enable_nudge::Bool
    nudge_time::TF
end


# Constructor for multidomain struct.
function MultiDomain(settings, TF)
    enable_nudge = settings["enable_nudge"]
    nudge_time = 1e12
    if enable_nudge
        nudge_time = settings["nudge_time"]
    end

    multidomain = MultiDomain{TF}(enable_nudge, nudge_time)
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


function calc_nudge_tend!(f_d::Fields, f_s::Fields, g_d::Grid, g_s::Grid, md::MultiDomain)
    if !md.enable_nudge
        return
    end

    c_nudge = 1 / md.nudge_time

    interp = LinearInterpolation((g_s.xh, g_s.y, g_s.z), f_s.u)
    f_d.u_tend[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] -= c_nudge * (f_d.u[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] - interp(g_d.xh[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke]))

    interp = LinearInterpolation((g_s.x, g_s.yh, g_s.z), f_s.v)
    f_d.v_tend[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] -= c_nudge * (f_d.v[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] - interp(g_d.x[g_d.is:g_d.ie], g_d.yh[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke]))

    interp = LinearInterpolation((g_s.x, g_s.y, g_s.zh), f_s.w)
    f_d.w_tend[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.keh] -= c_nudge * (f_d.w[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.keh] - interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.zh[g_d.ks:g_d.keh]))

    interp = LinearInterpolation((g_s.x, g_s.y, g_s.z), f_s.s)
    f_d.s_tend[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] -= c_nudge * (f_d.s[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] - interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke]))
end
