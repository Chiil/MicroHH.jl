using LoopVectorization

function advection_diffusion_kernel!(
    ut, u, v, w,
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

function advection_diffusion!(f::Fields, g::Grid)
    advection_diffusion_kernel!(
        f.u_tend, f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
