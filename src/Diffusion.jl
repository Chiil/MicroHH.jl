using .StencilBuilder


## Kernels.
function diffusion_u_kernel!(
    ut, u,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (ut[ih, j, k], u[ih, j, k]) begin
            ut += (
                + visc * (grad2x(grad2x(u)))
                + visc * (grad2y(grad2y(u)))
                + visc * (grad2z(grad2z(u))) )
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
                + visc * (grad2x(grad2x(v)))
                + visc * (grad2y(grad2y(v)))
                + visc * (grad2z(grad2z(v))) )
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
                + visc * (grad2x(grad2x(w)))
                + visc * (grad2y(grad2y(w)))
                + visc * (grad2z(grad2z(w))) )
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
                + visc * (grad2x(grad2x(s)))
                + visc * (grad2y(grad2y(s)))
                + visc * (grad2z(grad2z(s))) )
        end
    end
end


## Diffusion kernels launcher.
function calc_diffusion_tend!(f::Fields, g::Grid, to::TimerOutput)
    @timeit to "diffusion_u_kernel" diffusion_u_kernel!(
        f.u_tend, f.u,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    @timeit to "diffusion_v_kernel" diffusion_v_kernel!(
        f.v_tend, f.v,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    @timeit to "diffusion_w_kernel" diffusion_w_kernel!(
        f.w_tend, f.w,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        # g.is, g.ie, g.js, g.je, g.ks, g.keh) # CvH temporary hack
        g.is, g.ie, g.js, g.je, g.ks+1, g.keh-1)

    @timeit to "diffusion_s_kernel" diffusion_s_kernel!(
        f.s_tend, f.s,
        f.visc,
        g.dxi, g.dyi, g.dzi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    for scalar_name in keys(f.scalars)
        @timeit to "diffusion_s_kernel" diffusion_s_kernel!(
            f.scalars_tend[scalar_name], f.scalars[scalar_name],
            f.visc,
            g.dxi, g.dyi, g.dzi, g.dzhi,
            g.is, g.ie, g.js, g.je, g.ks, g.ke)
    end
end
