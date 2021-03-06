## Packages
include("../src/StencilBuilder.jl")

using BenchmarkTools
using LoopVectorization
using CUDA
using .StencilBuilder


## CPU dynamics kernel.
function dynamics_w_kernel_cpu!(
    wt, u, v, w, s,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    alpha = 9.81/300;

    @fast3d begin
        @fd (wt, u, v, w, s) wt += (
            - gradx(interpz(u) * interpx(w)) + visc * (gradx(gradx(w)))
            - grady(interpz(v) * interpy(w)) + visc * (grady(grady(w)))
            - gradz(interpz(w) * interpz(w)) + visc * (gradz(gradz(w)))
            + alpha*interpz(s) )
    end
end


## GPU dynamics kernel.
function dynamics_w_kernel_gpu!(
    wt, u, v, w, s,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    i = (blockIdx.().x - 1) * blockDim().x + threadIdx().x + is - 1
    j = (blockIdx.().y - 1) * blockDim().y + threadIdx().y + js - 1
    k = blockIdx.().z + ks - 1

    alpha = 9.81/300;

    if i <= ie && j <= je && k <= ke
        @inbounds @fd (wt, u, v, w, s) wt += (
            - gradx(interpz(u) * interpx(w)) + visc * (gradx(gradx(w)))
            - grady(interpz(v) * interpy(w)) + visc * (grady(grady(w)))
            - gradz(interpz(w) * interpz(w)) + visc * (gradz(gradz(w)))
            + alpha*interpz(s) )
    end

    return
end


## User input
itot = 384; jtot = 384; ktot = 384;
igc = 1; jgc = 1; kgc = 1;

icells = itot + 2igc; jcells = jtot + 2jgc; kcells = ktot + 2kgc
is = igc+1; ie = igc+itot; js = jgc+1; je = jgc+jtot; ks = kgc+1; ke = kgc+ktot

wt = rand(icells, jcells, kcells)
u = rand(icells, jcells, kcells)
v = rand(icells, jcells, kcells)
w = rand(icells, jcells, kcells)
s = rand(icells, jcells, kcells)
dxi = rand()
dyi = rand()
dzi = rand(kcells)
dzhi = rand(kcells)
visc = 1.


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
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    dynamics_w_kernel_cpu!(
        wt, u, v, w, s,
        visc,
        dxi, dyi, dzi, dzhi,
        is, ie, js, je, ks, ke)
end

function benchmark_gpu(
    wt_gpu, u_gpu, v_gpu, w_gpu, s_gpu,
    visc,
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
            visc,
            dxi, dyi, dzi_gpu, dzhi_gpu,
            is, ie, js, je, ks, ke)
    end
end


## Correctness test
benchmark_cpu(
    wt, u, v, w, s,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

benchmark_gpu(
    wt_gpu, u_gpu, v_gpu, w_gpu, s_gpu,
    visc,
    dxi, dyi, dzi_gpu, dzhi_gpu,
    is, ie, js, je, ks, ke)

wt_gpu_cpu = Array(wt_gpu)
println("Are wt and wt_gpu equal? ", isequal(wt, wt_gpu_cpu), " ", isapprox(wt, wt_gpu_cpu))


## Benchmarks
@btime benchmark_cpu(
    wt, u, v, w, s,
    visc,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

@btime benchmark_gpu(
    wt_gpu, u_gpu, v_gpu, w_gpu, s_gpu,
    visc,
    dxi, dyi, dzi_gpu, dzhi_gpu,
    is, ie, js, je, ks, ke)
