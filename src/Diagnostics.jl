function divergence_kernel(
    u, v, w,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    divmax = 0
    @fast3d begin
        div = @fd (u, v, w) gradx(u) + grady(v) + gradz(w)
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
        cfl = @fd (u, v, w) abs(interpx(u))*dxi + abs(interpy(v))*dyi + abs(interpz(w))*dzi[k]
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
