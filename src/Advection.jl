using .StencilBuilder


## Kernels.
function advection_u_kernel!(
    ut, u, v, w,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (ut[ih, j, k], u[ih, j, k], v[i, jh, k], w[i, j, kh]) begin
            ut += (
                - grad2x(interp2x(u) * interp2x(u))
                - grad2y(interp2x(v) * interp2y(u))
                - grad2z(interp2x(w) * interp2z(u)) )
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
                - grad2x(interp2y(u) * interp2x(v))
                - grad2y(interp2y(v) * interp2y(v))
                - grad2z(interp2y(w) * interp2z(v)) )
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
                - grad2x(interp2z(u) * interp2x(w))
                - grad2y(interp2z(v) * interp2y(w))
                - grad2z(interp2z(w) * interp2z(w)) )
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
                - grad2x(u * interp2x(s))
                - grad2y(v * interp2y(s))
                - grad2z(w * interp2z(s)) )
        end
    end
end


## Advection kernels launcher.
function calc_advection_tend!(f::Fields, g::Grid, to::TimerOutput)
    @timeit to "advection_u_kernel" advection_u_kernel!(
        f.u_tend, f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    @timeit to "advection_v_kernel" advection_v_kernel!(
        f.v_tend, f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    @timeit to "advection_w_kernel" advection_w_kernel!(
        f.w_tend, f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzhi,
        # g.is, g.ie, g.js, g.je, g.ks, g.keh) # CvH temporary hack
        g.is, g.ie, g.js, g.je, g.ks+1, g.keh-1)

    @timeit to "advection_s_kernel" advection_s_kernel!(
        f.s_tend, f.u, f.v, f.w, f.s,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    for scalar_name in keys(f.scalars)
        @timeit to "advection_s_kernel" advection_s_kernel!(
            f.scalars_tend[scalar_name], f.u, f.v, f.w, f.scalars[scalar_name],
            g.dxi, g.dyi, g.dzi,
            g.is, g.ie, g.js, g.je, g.ks, g.ke)
    end
end
