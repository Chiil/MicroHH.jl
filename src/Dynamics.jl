using .StencilBuilder


## Kernels.
function advection_u_kernel!(
    ut, u, v, w,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (ut[ih, j, k], u[ih, j, k], v[i, jh, k], w[i, j, kh]) begin
            ut += (
                - gradx(interpx(u) * interpx(u))
                - grady(interpx(v) * interpy(u))
                - gradz(interpx(w) * interpz(u)) )
        end
    end
end


function diffusion_u_kernel!(
    ut, u,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (ut[ih, j, k], u[ih, j, k]) begin
            ut += (
                + visc * (gradx(gradx(u)))
                + visc * (grady(grady(u)))
                + visc * (gradz(gradz(u))) )
        end
    end
end


function advection_v_kernel!(
    vt, u, v, w,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (vt[i, jh, k], u[ih, j, k], v[i, jh, k], w[i, j, kh]) begin
            vt += (
                - gradx(interpy(u) * interpx(v))
                - grady(interpy(v) * interpy(v))
                - gradz(interpy(w) * interpz(v)) )
        end
    end
end


function diffusion_v_kernel!(
    vt, v,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (vt[i, jh, k], v[i, jh, k]) begin
            vt += (
                + visc * (gradx(gradx(v)))
                + visc * (grady(grady(v)))
                + visc * (gradz(gradz(v))) )
        end
    end
end


function advection_w_kernel!(
    wt, u, v, w,
    dxi, dyi, dzhi,
    is, ie, js, je, ks, ke) # CvH: ke represents keh here, but that is not yet supported by fast3d, needs fixing.

    @fast3d begin
        @fd (wt[i, j, kh], u[ih, j, k], v[i, jh, k], w[i, j, kh]) begin
            wt += (
                - gradx(interpz(u) * interpx(w))
                - grady(interpz(v) * interpy(w))
                - gradz(interpz(w) * interpz(w)) )
        end
    end
end


function diffusion_w_kernel!(
    wt, w,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke) # CvH: ke represents keh here, but that is not yet supported by fast3d, needs fixing.

    @fast3d begin
        @fd (wt[i, j, kh], w[i, j, kh]) begin
            wt += (
                + visc * (gradx(gradx(w)))
                + visc * (grady(grady(w)))
                + visc * (gradz(gradz(w))) )
        end
    end
end


function buoyancy_w_kernel!(
    wt, s,
    s_ref, alpha,
    is, ie, js, je, ks, ke) # CvH: ke represents keh here, but that is not yet supported by fast3d, needs fixing.

    @fast3d begin
        @fd (wt[i, j, kh], s[i, j, k], s_ref[k]) begin
            wt += alpha*interpz(s - s_ref)
        end
    end
end


function advection_s_kernel!(
    st, u, v, w, s,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (st[i, j, k], u[ih, j, k], v[i, jh, k], w[i, j, kh], s[i, j, k]) begin
            st += (
                - gradx(u * interpx(s))
                - grady(v * interpy(s))
                - gradz(w * interpz(s)) )
        end
    end
end


function diffusion_s_kernel!(
    st, s,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (st[i, j, k], s[i, j, k]) begin
            st += (
                + visc * (gradx(gradx(s)))
                + visc * (grady(grady(s)))
                + visc * (gradz(gradz(s))) )
        end
    end
end


## Dynamics kernels launcher.
function calc_dynamics_tend!(f::Fields, g::Grid, to::TimerOutput)
    @timeit to "advection_u_kernel" advection_u_kernel!(
        f.u_tend, f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    @timeit to "diffusion_u_kernel" diffusion_u_kernel!(
        f.u_tend, f.u,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    @timeit to "advection_v_kernel" advection_v_kernel!(
        f.v_tend, f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    @timeit to "diffusion_v_kernel" diffusion_v_kernel!(
        f.v_tend, f.v,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    @timeit to "advection_w_kernel" advection_w_kernel!(
        f.w_tend, f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzhi,
        # g.is, g.ie, g.js, g.je, g.ks, g.keh) # CvH temporary hack
        g.is, g.ie, g.js, g.je, g.ks+1, g.keh-1)

    @timeit to "diffusion_w_kernel" diffusion_w_kernel!(
        f.w_tend, f.w,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        # g.is, g.ie, g.js, g.je, g.ks, g.keh) # CvH temporary hack
        g.is, g.ie, g.js, g.je, g.ks+1, g.keh-1)

    @timeit to "buoyancy_w_kernel" buoyancy_w_kernel!(
        f.w_tend, f.s,
        f.s_ref, f.alpha,
        # g.is, g.ie, g.js, g.je, g.ks, g.keh) # CvH temporary hack
        g.is, g.ie, g.js, g.je, g.ks+1, g.keh-1)

    @timeit to "advection_s_kernel" advection_s_kernel!(
        f.s_tend, f.u, f.v, f.w, f.s,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    @timeit to "diffusion_s_kernel" diffusion_s_kernel!(
        f.s_tend, f.s,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
