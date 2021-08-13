function input_kernel!(
    p,
    u, v, w,
    ut, vt, wt,
    dxi, dyi, dzi, dti,
    is, ie, js, je, ks, ke)

    @fast3d begin
        p[i-is+1, j-js+1, k-ks+1] = @fd (u, v, w, ut, vt, wt) gradx(ut + u*dti) + grady(vt + v*dti) + gradz(wt + w*dti)
    end
end

function calc_pressure_tend!(f::Fields, g::Grid, t::Timeloop)
    input_kernel!(
        f.p,
        f.u, f.v, f.w,
        f.u_tend, f.v_tend, f.w_tend,
        g.dxi, g.dyi, g.dzi, 1/t.dt,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
