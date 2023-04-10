using .StencilBuilder


## Kernels.
function buoyancy_w_kernel!(
    wt, s,
    s_ref, alpha,
    is, ie, js, je, ks, ke) # CvH: ke represents keh here, but that is not yet supported by fast3d, needs fixing.

    @fast3d begin
        @fd (wt[i, j, kh], s[i, j, k], s_ref[k]) begin
            wt += alpha*interp2z(s - s_ref)
        end
    end
end


## Buoyancy kernel launcher.
function calc_buoyancy_tend!(f::Fields, g::Grid, to::TimerOutput)
    @timeit to "buoyancy_w_kernel" buoyancy_w_kernel!(
        f.w_tend, f.s,
        f.s_ref, f.alpha,
        # g.is, g.ie, g.js, g.je, g.ks, g.keh) # CvH temporary hack
        g.is, g.ie, g.js, g.je, g.ks+1, g.keh-1)
end
