using .StencilBuilder


## Kernels.
function dynamics_u_kernel!(
    ut, u, v, w,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (ut, u, v, w) ut += (
            - gradx(interpx(u) * interpx(u)) + visc * (gradx(gradx(u)))
            - grady(interpx(v) * interpy(u)) + visc * (grady(grady(u)))
            - gradz(interpx(w) * interpz(u)) + visc * (gradz(gradz(u))) )
    end
end


function dynamics_v_kernel!(
    vt, u, v, w,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (vt, u, v, w) vt += (
            - gradx(interpy(u) * interpx(v)) + visc * (gradx(gradx(v)))
            - grady(interpy(v) * interpy(v)) + visc * (grady(grady(v)))
            - gradz(interpy(w) * interpz(v)) + visc * (gradz(gradz(v))) )
    end
end


function dynamics_w_kernel!(
    wt, u, v, w, s,
    visc, alpha,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (wt, u, v, w, s) wt += (
            - gradx(interpz(u) * interpx(w)) + visc * (gradx(gradx(w)))
            - grady(interpz(v) * interpy(w)) + visc * (grady(grady(w)))
            - gradz(interpz(w) * interpz(w)) + visc * (gradz(gradz(w)))
            + alpha*interpz(s) )
    end
end


function dynamics_s_kernel!(
    st, u, v, w, s,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (st, u, v, w, s) st += (
            - gradx(u * interpx(s)) + visc * (gradx(gradx(s)))
            - grady(v * interpy(s)) + visc * (grady(grady(s)))
            - gradz(w * interpz(s)) + visc * (gradz(gradz(s))) )
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
        f.w_tend, f.u, f.v, f.w, f.s,
        f.visc, f.alpha,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    dynamics_s_kernel!(
        f.s_tend, f.u, f.v, f.w, f.s,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
