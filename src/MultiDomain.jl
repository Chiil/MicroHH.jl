using Interpolations

function transfer_state(f_d, f_s, g_d, g_s)
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