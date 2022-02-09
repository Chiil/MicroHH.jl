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
    interp = interpolate((g_s.xh, g_s.y, g_s.z), f_s.u, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    f_d.u[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] .= interp(g_d.xh[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = interpolate((g_s.x, g_s.yh, g_s.z), f_s.v, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    f_d.v[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] .= interp(g_d.x[g_d.is:g_d.ie], g_d.yh[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = interpolate((g_s.x, g_s.y, g_s.zh), f_s.w, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    f_d.w[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.keh] .= interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.zh[g_d.ks:g_d.keh])

    interp = interpolate((g_s.x, g_s.y, g_s.z), f_s.s, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    f_d.s[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] .= interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = interpolate((g_s.x, g_s.y, g_s.z), f_s.p, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    f_d.p[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] .= interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])

    interp = interpolate((g_s.x, g_s.y), f_s.s_bot, (Gridded(Linear()), Gridded(Linear())))
    f_d.s_bot[g_d.is:g_d.ie, g_d.js:g_d.je] .= interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je])

    interp = interpolate((g_s.x, g_s.y), f_s.s_gradbot, (Gridded(Linear()), Gridded(Linear())))
    f_d.s_gradbot[g_d.is:g_d.ie, g_d.js:g_d.je] .= interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je])

    interp = interpolate((g_s.x, g_s.y), f_s.s_top, (Gridded(Linear()), Gridded(Linear())))
    f_d.s_top[g_d.is:g_d.ie, g_d.js:g_d.je] .= interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je])

    interp = interpolate((g_s.x, g_s.y), f_s.s_gradtop, (Gridded(Linear()), Gridded(Linear())))
    f_d.s_gradtop[g_d.is:g_d.ie, g_d.js:g_d.je] .= interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je])
end


## Define the upsample functions.
function make_inner_loop!(ex_list, ii, jj, kk, ifac, jfac, kfac, ioff, joff, koff, is_top)
    push!(ex_list, :(i_hi = $ifac*i+$ii+is_hi), :(j_hi = $jfac*j+$jj+js_hi))
    if is_top
        push!(ex_list, :(k_hi = ke_hi))
    else
        push!(ex_list, :(k_hi = $kfac*k+$kk+ks_hi))
    end

    i_pos = -1/2 + ioff + (1/2 - ioff + ii)/ifac
    j_pos = -1/2 + joff + (1/2 - joff + jj)/jfac
    k_pos = -1/2 + koff + (1/2 - koff + kk)/kfac

    i_lo_off = floor(Int, i_pos)
    j_lo_off = floor(Int, j_pos)
    k_lo_off = floor(Int, k_pos)

    if is_top
        push!(ex_list, :(i_lo = i + $i_lo_off + is_lo), :(j_lo = j + $j_lo_off + js_lo), :(k_lo = ke_lo))
    else
        push!(ex_list, :(i_lo = i + $i_lo_off + is_lo), :(j_lo = j + $j_lo_off + js_lo), :(k_lo = k + $k_lo_off + ks_lo))
    end

    fi = mod(i_pos, 1)
    fj = mod(j_pos, 1)
    fk = mod(k_pos, 1)

    coef = zeros(2, 2, 2)
    coef[1, 1, 1] = (((1-fi)*fj+fi-1)*fk+(fi-1)*fj-fi+1)
    coef[2, 1, 1] = ((fi*fj-fi)*fk-fi*fj+fi)
    coef[1, 2, 1] = ((fi-1)*fj*fk+(1-fi)*fj)
    coef[2, 2, 1] = (fi*fj-fi*fj*fk)
    coef[1, 1, 2] = ((fi-1)*fj-fi+1)*fk
    coef[2, 1, 2] = (fi-fi*fj)*fk
    coef[1, 2, 2] = (1-fi)*fj*fk
    coef[2, 2, 2] = fi*fj*fk

    kc_end = is_top ? 1 : 2

    rhs_list = []
    for kc in 1:kc_end, jc in 1:2, ic in 1:2
        if !(coef[ic, jc, kc] ≈ 0)
            push!(rhs_list, :( $(coef[ic, jc, kc])*lo[i_lo+$(ic-1), j_lo+$(jc-1), k_lo+$(kc-1)] ))
        end
    end

    rhs = Expr(:call, :+, rhs_list...)

    push!(ex_list, :(hi[i_hi, j_hi, k_hi] = $rhs))
end


macro upsample_lin_fast(suffix, ifac, jfac, kfac, ioff, joff, koff, add_top)
    ex_inner = []
    for kk in 0:kfac-1, jj in 0:jfac-1, ii in 0:ifac-1
        make_inner_loop!(ex_inner, ii, jj, kk, ifac, jfac, kfac, ioff, joff, koff, false)
    end
    ex_inner_block = quote
        Threads.@threads for k in 0:ktot-1
            for j in 0:jmax-1
                @inbounds @fastmath @simd for i in 0:imax-1
                    $(Expr(:block, ex_inner...))
                end
            end
        end
    end

    if add_top
        ex_top = []
        kk = 0
        for jj in 0:jfac-1, ii in 0:ifac-1
            make_inner_loop!(ex_top, ii, jj, kk, ifac, jfac, kfac, ioff, joff, koff, true)
        end
        ex_top_block = quote
            Threads.@threads for j in 0:jmax-1
                @inbounds @fastmath @simd for i in 0:imax-1
                    $(Expr(:block, ex_top...))
                end
            end
        end
    else
        ex_top_block = :()
    end

    name = Symbol(@sprintf "upsample_lin_%s!" suffix)
    ex = quote
        function $name(hi, lo, imax, jmax, ktot, is_hi, js_hi, ks_hi, ke_hi, is_lo, js_lo, ks_lo, ke_lo)
            $ex_inner_block
            $ex_top_block
        end
    end

    return esc(ex)
end


@upsample_lin_fast("222_u", 2, 2, 2, 0.5, 0  ,   0, false)
@upsample_lin_fast("222_v", 2, 2, 2,   0, 0.5,   0, false)
@upsample_lin_fast("222_w", 2, 2, 2,   0, 0  , 0.5, true )
@upsample_lin_fast("222_s", 2, 2, 2,   0, 0  ,   0, false)
@upsample_lin_fast("333_u", 3, 3, 3, 0.5, 0  ,   0, false)
@upsample_lin_fast("333_v", 3, 3, 3,   0, 0.5,   0, false)
@upsample_lin_fast("333_w", 3, 3, 3,   0, 0  , 0.5, true )
@upsample_lin_fast("333_s", 3, 3, 3,   0, 0  ,   0, false)
@upsample_lin_fast("444_u", 4, 4, 4, 0.5, 0  ,   0, false)
@upsample_lin_fast("444_v", 4, 4, 4,   0, 0.5,   0, false)
@upsample_lin_fast("444_w", 4, 4, 4,   0, 0  , 0.5, true )
@upsample_lin_fast("444_s", 4, 4, 4,   0, 0  ,   0, false)


function calc_nudge_fields!(md::MultiDomain, f_d::Fields, f_s::Fields, g_d::Grid, g_s::Grid)
    if !md.enable_nudge
        return
    end

    upsample_ratio = (g_d.itot, g_d.jtot, g_d.ktot) .÷ (g_s.itot, g_s.jtot, g_s.ktot)
    if upsample_ratio == (2, 2, 2)
        upsample_lin_222_u!(
            md.u_nudge, f_s.u, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

        upsample_lin_222_v!(
            md.v_nudge, f_s.v, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

        upsample_lin_222_w!(
            md.w_nudge, f_s.w, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.keh, g_s.is, g_s.js, g_s.ks, g_s.keh)

        upsample_lin_222_s!(
            md.s_nudge, f_s.s, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

    elseif upsample_ratio == (3, 3, 3)
        upsample_lin_333_u!(
            md.u_nudge, f_s.u, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

        upsample_lin_333_v!(
            md.v_nudge, f_s.v, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

        upsample_lin_333_w!(
            md.w_nudge, f_s.w, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.keh, g_s.is, g_s.js, g_s.ks, g_s.keh)

        upsample_lin_333_s!(
            md.s_nudge, f_s.s, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

    elseif upsample_ratio == (4, 4, 4)
        upsample_lin_444_u!(
            md.u_nudge, f_s.u, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

        upsample_lin_444_v!(
            md.v_nudge, f_s.v, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

        upsample_lin_444_w!(
            md.w_nudge, f_s.w, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.keh, g_s.is, g_s.js, g_s.ks, g_s.keh)

        upsample_lin_444_s!(
            md.s_nudge, f_s.s, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

    elseif upsample_ratio == (4, 4, 4)
        upsample_lin_444_u!(
            md.u_nudge, f_s.u, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

        upsample_lin_444_v!(
            md.v_nudge, f_s.v, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

        upsample_lin_444_w!(
            md.w_nudge, f_s.w, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.keh, g_s.is, g_s.js, g_s.ks, g_s.keh)

        upsample_lin_444_s!(
            md.s_nudge, f_s.s, g_s.imax, g_s.jmax, g_s.ktot,
            g_d.is, g_d.js, g_d.ks, g_d.ke, g_s.is, g_s.js, g_s.ks, g_s.ke)

    else
        println("WARNING: resorting to slow interpolations in nudging.")
        @sync begin
            Threads.@spawn begin
                interp = interpolate((g_s.xh, g_s.y, g_s.z), f_s.u, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
                md.u_nudge[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] .= interp(g_d.xh[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])
            end

            Threads.@spawn begin
                interp = interpolate((g_s.x, g_s.yh, g_s.z), f_s.v, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
                md.v_nudge[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] .= interp(g_d.x[g_d.is:g_d.ie], g_d.yh[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])
            end

            Threads.@spawn begin
                interp = interpolate((g_s.x, g_s.y, g_s.zh), f_s.w, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
                md.w_nudge[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.keh] .= interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.zh[g_d.ks:g_d.keh])
            end

            Threads.@spawn begin
                interp = interpolate((g_s.x, g_s.y, g_s.z), f_s.s, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
                md.s_nudge[g_d.is:g_d.ie, g_d.js:g_d.je, g_d.ks:g_d.ke] .= interp(g_d.x[g_d.is:g_d.ie], g_d.y[g_d.js:g_d.je], g_d.z[g_d.ks:g_d.ke])
            end
        end
    end
end


function nudge_kernel!(a_tend, a, a_nudge, c_nudge, is, ie, js, je, ks, ke)
    @fast3d begin
        @fd (a_tend[i, j, k], a[i, j, k], a_nudge[i, j, k]) a_tend -= c_nudge * (a - a_nudge)
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
