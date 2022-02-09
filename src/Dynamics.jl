using .StencilBuilder


## Kernels.
function dynamics_u_kernel!(
    ut, u, v, w,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (ut[ih, j, k], u[ih, j, k], v[i, jh, k], w[i, j, kh]) begin
            ut += (
                - gradx(interpx(u) * interpx(u)) + visc * (gradx(gradx(u)))
                - grady(interpx(v) * interpy(u)) + visc * (grady(grady(u)))
                - gradz(interpx(w) * interpz(u)) + visc * (gradz(gradz(u))) )
        end
    end
end


function dynamics_v_kernel!(
    vt, u, v, w,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (vt[i, jh, k], u[ih, j, k], v[i, jh, k], w[i, j, kh]) begin
            vt += (
                - gradx(interpy(u) * interpx(v)) + visc * (gradx(gradx(v)))
                - grady(interpy(v) * interpy(v)) + visc * (grady(grady(v)))
                - gradz(interpy(w) * interpz(v)) + visc * (gradz(gradz(v))) )
        end
    end
end


function dynamics_w_kernel!(
    wt, u, v, w, s, s_ref,
    visc, alpha,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke) # CvH: ke represents keh here, but that is not yet supported by fast3d, needs fixing.

    @fast3d begin
        @fd (wt[i, j, kh], u[ih, j, k], v[i, jh, k], w[i, j, kh], s[i, j, k], s_ref[k]) begin
            wt += (
                - gradx(interpz(u) * interpx(w)) + visc * (gradx(gradx(w)))
                - grady(interpz(v) * interpy(w)) + visc * (grady(grady(w)))
                - gradz(interpz(w) * interpz(w)) + visc * (gradz(gradz(w)))
                + alpha*interpz(s - s_ref) )
        end
    end
end


function dynamics_s_kernel!(
    st, u, v, w, s,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (st[i, j, k], u[ih, j, k], v[i, jh, k], w[i, j, kh], s[i, j, k]) begin
            st += (
                - gradx(u * interpx(s)) + visc * (gradx(gradx(s)))
                - grady(v * interpy(s)) + visc * (grady(grady(s)))
                - gradz(w * interpz(s)) + visc * (gradz(gradz(s))) )
        end
    end
end


## Dynamics kernels launcher.
function calc_dynamics_tend!(f::Fields, g::Grid)
    dynamics_u_kernel!(
        f.u_tend, f.u, f.v, f.w,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    dynamics_v_kernel!(
        f.v_tend, f.u, f.v, f.w,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    dynamics_w_kernel!(
        f.w_tend, f.u, f.v, f.w, f.s, f.s_ref,
        f.visc, f.alpha,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        # g.is, g.ie, g.js, g.je, g.ks, g.keh) # CvH temporary hack
        g.is, g.ie, g.js, g.je, g.ks+1, g.keh-1)

    dynamics_s_kernel!(
        f.s_tend, f.u, f.v, f.w, f.s,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
