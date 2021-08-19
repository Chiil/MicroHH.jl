## Packages
include("../src/StencilBuilder.jl")

using BenchmarkTools
using LoopVectorization
using .StencilBuilder


## CPU dynamics kernel.
function dynamics_w_kernel!(
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


## Run CPU kernel
@btime dynamics_w_kernel!(
    $wt, $u, $v, $w, $s,
    $visc,
    $dxi, $dyi, $dzi, $dzhi,
    $is, $ie, $js, $je, $ks, $ke)

