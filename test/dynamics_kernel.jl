## Packages and settings.
using BenchmarkTools
using LoopVectorization
using MicroHH.StencilBuilder

float_type = Float32


## Dynamics kernel with all processes combined.
function dynamics_w_kernel!(
    wt, u, v, w, s,
    visc, alpha,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (wt[i, j, kh], u[ih, j, k], v[i, jh, k], w[i, j, kh], s[i, j, k]) begin
            wt += (
                - grad2x(interp2z(u) * interp2x(w)) + visc * (grad2x(grad2x(w)))
                - grad2y(interp2z(v) * interp2y(w)) + visc * (grad2y(grad2y(w)))
                - grad2z(interp2z(w) * interp2z(w)) + visc * (grad2z(grad2z(w)))
                + alpha*interp2z(s) )
        end
    end
end


## Dynamics kernel split per process.
function dynamics_w_kernel_1!(
    wt, u, v, w,
    dxi, dyi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (wt[i, j, kh], u[ih, j, k], v[i, jh, k], w[i, j, kh]) begin
            wt += (
                - grad2x(interp2z(u) * interp2x(w))
                - grad2y(interp2z(v) * interp2y(w))
                - grad2z(interp2z(w) * interp2z(w)) )
        end
    end
end


function dynamics_w_kernel_2!(
    wt, w,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (wt[i, j, kh], w[i, j, kh]) begin
            wt += (
                + visc * (grad2x(grad2x(w)))
                + visc * (grad2y(grad2y(w)))
                + visc * (grad2z(grad2z(w))) )
        end
    end
end


function dynamics_w_kernel_3!(
    wt, s,
    alpha,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (wt[i, j, kh], s[i, j, k]) begin
            wt += alpha*interp2z(s)
        end
    end
end


## User input
itot = 384; jtot = 384; ktot = 384;
igc = 1; jgc = 1; kgc = 1;

icells = itot + 2igc; jcells = jtot + 2jgc; kcells = ktot + 2kgc
is = igc+1; ie = igc+itot; js = jgc+1; je = jgc+jtot; ks = kgc+1; ke = kgc+ktot

wt = rand(float_type, icells, jcells, kcells)
u = rand(float_type, icells, jcells, kcells)
v = rand(float_type, icells, jcells, kcells)
w = rand(float_type, icells, jcells, kcells)
s = rand(float_type, icells, jcells, kcells)
dxi = rand(float_type)
dyi = rand(float_type)
dzi = rand(float_type, kcells)
dzhi = rand(float_type, kcells)
visc = convert(float_type, 1)
alpha = convert(float_type, 9.81/300)


## Benchmark kernels.
@btime dynamics_w_kernel!(
    $wt, $u, $v, $w, $s,
    $visc, $alpha,
    $dxi, $dyi, $dzi, $dzhi,
    $is, $ie, $js, $je, $ks, $ke)

@btime begin
    dynamics_w_kernel_1!(
        $wt, $u, $v, $w,
        $dxi, $dyi, $dzhi,
        $is, $ie, $js, $je, $ks, $ke)

    dynamics_w_kernel_2!(
        $wt, $w,
        $visc,
        $dxi, $dyi, $dzi, $dzhi,
        $is, $ie, $js, $je, $ks, $ke)

    dynamics_w_kernel_3!(
        $wt, $s,
        $alpha,
        $is, $ie, $js, $je, $ks, $ke)
end

@btime begin
    dynamics_w_kernel_1!(
        $wt, $u, $v, $w,
        $dxi, $dyi, $dzhi,
        $is, $ie, $js, $je, $ks, $ke)
end

@btime begin
    dynamics_w_kernel_2!(
        $wt, $w,
        $visc,
        $dxi, $dyi, $dzi, $dzhi,
        $is, $ie, $js, $je, $ks, $ke)
end

@btime begin
    dynamics_w_kernel_3!(
        $wt, $s,
        $alpha,
        $is, $ie, $js, $je, $ks, $ke)
end
