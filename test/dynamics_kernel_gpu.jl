## Packages
include("../src/StencilBuilder.jl")

using BenchmarkTools
using LoopVectorization
using Tullio
using CUDA
using .StencilBuilder


## CPU dynamics kernel.
function dynamics_w_kernel_cpu!(
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


## GPU dynamics kernel.
macro calc_ijk(is, js, ks)
    i = :( (blockIdx.().x - 1) * blockDim().x + threadIdx().x + $is - 1 )
    j = :( (blockIdx.().y - 1) * blockDim().y + threadIdx().y + $js - 1 )
    k = :( blockIdx.().z + $ks - 1 )

    ex = :(i = $i; j = $j; k = $k)
    return(esc(ex))
end

function dynamics_w_kernel_gpu!(
    wt, u, v, w, s,
    visc, alpha,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @calc_ijk(is, js, ks)

    if i <= ie && j <= je && k <= ke
        @inbounds @fd (wt[i, j, kh], u[ih, j, k], v[i, jh, k], w[i, j, kh], s[i, j, k]) begin
            wt += (
                - grad2x(interp2z(u) * interp2x(w)) + visc * (grad2x(grad2x(w)))
                - grad2y(interp2z(v) * interp2y(w)) + visc * (grad2y(grad2y(w)))
                - grad2z(interp2z(w) * interp2z(w)) + visc * (grad2z(grad2z(w)))
                + alpha*interp2z(s) )
        end
    end

    return
end


## User input
float_type = Float64
itot = 768; jtot = 768; ktot = 384;
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
visc = 1.
alpha = 9.81/300.


## Copy data to GPU.
wt_gpu = CuArray(wt)
u_gpu = CuArray(u)
v_gpu = CuArray(v)
w_gpu = CuArray(w)
s_gpu = CuArray(s)
dzi_gpu = CuArray(dzi)
dzhi_gpu = CuArray(dzhi)


## Benchmark functions
function benchmark_cpu(
    wt, u, v, w, s,
    visc, alpha,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    dynamics_w_kernel_cpu!(
        wt, u, v, w, s,
        visc, alpha,
        dxi, dyi, dzi, dzhi,
        is, ie, js, je, ks, ke)
end

function benchmark_gpu(
    wt_gpu, u_gpu, v_gpu, w_gpu, s_gpu,
    visc, alpha,
    dxi, dyi, dzi_gpu, dzhi_gpu,
    is, ie, js, je, ks, ke)

    blocks_size = (32, 32, 1)
    blocks_num = (
        ceil(Int, (ie-is+1) / blocks_size[1]),
        ceil(Int, (je-js+1) / blocks_size[2]),
        (ke-ks+1) )

    CUDA.@sync begin
        @cuda threads=blocks_size blocks=blocks_num dynamics_w_kernel_gpu!(
            wt_gpu, u_gpu, v_gpu, w_gpu, s_gpu,
            visc, alpha,
            dxi, dyi, dzi_gpu, dzhi_gpu,
            is, ie, js, je, ks, ke)
    end
end

function benchmark_tullio(
    wt, u, v, w, s,
    visc, alpha,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fd_tullio (wt[i, j, kh], u[ih, j, k], v[i, jh, k], w[i, j, kh], s[i, j, k]) begin
        wt += (
            - grad2x(interp2z(u) * interp2x(w)) + visc * (grad2x(grad2x(w)))
            - grad2y(interp2z(v) * interp2y(w)) + visc * (grad2y(grad2y(w)))
            - grad2z(interp2z(w) * interp2z(w)) + visc * (grad2z(grad2z(w)))
            + alpha*interp2z(s) )
    end

end


## Correctness test
benchmark_cpu(
    wt, u, v, w, s,
    visc, alpha,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

# benchmark_gpu(
benchmark_gpu(
    wt_gpu, u_gpu, v_gpu, w_gpu, s_gpu,
    visc, alpha,
    dxi, dyi, dzi_gpu, dzhi_gpu,
    is, ie, js, je, ks, ke)

wt_gpu_cpu = Array(wt_gpu)
println("Are wt and wt_gpu (exactly, approximately) equal? ", isequal(wt, wt_gpu_cpu), " ", isapprox(wt, wt_gpu_cpu))


## Benchmarks
@btime benchmark_cpu(
    $wt, $u, $v, $w, $s,
    $visc, $alpha,
    $dxi, $dyi, $dzi, $dzhi,
    $is, $ie, $js, $je, $ks, $ke)

@btime benchmark_gpu(
    $wt_gpu, $u_gpu, $v_gpu, $w_gpu, $s_gpu,
    $visc, $alpha,
    $dxi, $dyi, $dzi_gpu, $dzhi_gpu,
    $is, $ie, $js, $je, $ks, $ke)

@btime benchmark_tullio(
    # wt_gpu, u_gpu, v_gpu, w_gpu, s_gpu,
    # visc, alpha,
    # dxi, dyi, dzi_gpu, dzhi_gpu,
    # is, ie, js, je, ks, ke)
    $wt, $u, $v, $w, $s,
    $visc, $alpha,
    $dxi, $dyi, $dzi, $dzhi,
    $is, $ie, $js, $je, $ks, $ke)
