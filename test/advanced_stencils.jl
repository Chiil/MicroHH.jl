## Packages.
using Tullio
using LoopVectorization

include("../src/StencilBuilder.jl")
using .StencilBuilder


## Settings.
itot = 128; jtot = 96; ktot = 64
igc = 1; jgc = 1; kgc = 2
is = igc+1; js = jgc+1; ks = kgc+1
ie = igc+itot; je = jgc+jtot; ke = kgc+ktot

alpha = 1.
u = rand(itot+2igc, jtot+2jgc, ktot+2kgc)
v = rand(itot+2igc, jtot+2jgc, ktot+2kgc)
w = rand(itot+2igc, jtot+2jgc, ktot+2kgc)
s = rand(itot+2igc, jtot+2jgc, ktot+2kgc)
s_ref = rand(ktot+2kgc)

wt = rand(itot+2igc, jtot+2jgc, ktot+2kgc)


## Settings.
# @fd_tullio () (0, 0, -1/2) tmp = alpha * interpz(s - s_ref)
# @fd_tullio () wt = alpha * interpz(s - s_ref)
@fast3d begin
    @fd (wt[i, j, kh], u[ih, j, k], v[i, jh, k], w[i, j, kh], s[i, j, k], s_ref[k]) begin
        wt += (
            - gradx(interpz(u) * interpx(w)) + visc * (gradx(gradx(w)))
            - grady(interpz(v) * interpy(w)) + visc * (grady(grady(w)))
            - gradz(interpz(w) * interpz(w)) + visc * (gradz(gradz(w)))
            + alpha*interpz(s - s_ref) )
    end
end
