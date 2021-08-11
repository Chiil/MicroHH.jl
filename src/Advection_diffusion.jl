using LoopVectorization


function advec_diff_u_kernel!(
    ut, u, v, w,
    visc,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @tturbo for k in ks:ke
        for j in js:je
            for i in is:ie
                ut[i, j, k] += u[i, j, k] + v[i, j, k] + w[i, j, k]
            end
        end
    end
end


function advec_diff_v_kernel!(
    vt, u, v, w,
    visc,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @tturbo for k in ks:ke
        for j in js:je
            for i in is:ie
                vt[i, j, k] += u[i, j, k] + v[i, j, k] + w[i, j, k]
            end
        end
    end
end


function advec_diff_w_kernel!(
    wt, u, v, w,
    visc,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @tturbo for k in ks:ke
        for j in js:je
            for i in is:ie
                wt[i, j, k] += u[i, j, k] + v[i, j, k] + w[i, j, k]
            end
        end
    end
end


function advec_diff_s_kernel!(
    st, u, v, w, s,
    visc,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @tturbo for k in ks:ke
        for j in js:je
            for i in is:ie
                st[i, j, k] += u[i, j, k] + v[i, j, k] + w[i, j, k]
            end
        end
    end
end


function advection_diffusion!(f::Fields, g::Grid)
    advec_diff_u_kernel!(
        f.u_tend, f.u, f.v, f.w,
        f.visc,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    advec_diff_v_kernel!(
        f.v_tend, f.u, f.v, f.w,
        f.visc,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    advec_diff_w_kernel!(
        f.w_tend, f.u, f.v, f.w,
        f.visc,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    advec_diff_s_kernel!(
        f.s_tend, f.u, f.v, f.w, f.s,
        f.visc,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
