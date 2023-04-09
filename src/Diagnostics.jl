function divergence_kernel(
    u, v, w,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    divmax = 0
    @fast3d begin
        @fd (u[ih, j, k], v[i, jh, k], w[i, j, kh]) begin
            div = grad2x(u) + grad2y(v) + grad2z(w)
        end
        divmax = max(divmax, abs(div))
    end

    return divmax
end


function cfl_kernel(
    u, v, w,
    dxi, dyi, dzi, dt,
    is, ie, js, je, ks, ke)

    cflmax = 0
    @fast3d begin
        @fd (u[ih, j, k], v[i, jh, k], w[i, j, kh], dzi[k]) begin
            cfl = abs(interp2x(u))*dxi + abs(interp2y(v))*dyi + abs(interp2z(w))*dzi
        end
        cflmax = max(cflmax, cfl)
    end

    return cflmax*dt
end


function calc_divergence(f::Fields, g::Grid)
    div = divergence_kernel(
        f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end


function calc_cfl(f::Fields, g::Grid, t::Timeloop)
    cfl = cfl_kernel(
        f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzi, t.dt,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
