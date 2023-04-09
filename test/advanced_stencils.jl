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
# @fd_tullio () (0, 0, -1/2) tmp = alpha * interp2z(s - s_ref)
# @fd_tullio () wt = alpha * interp2z(s - s_ref)
@fast3d begin
    @fd (wt[i, j, kh], u[ih, j, k], v[i, jh, k], w[i, j, kh], s[i, j, k], s_ref[k]) begin
        wt += (
            - grad2x(interp2z(u) * interp2x(w)) + visc * (grad2x(grad2x(w)))
            - grad2y(interp2z(v) * interp2y(w)) + visc * (grad2y(grad2y(w)))
            - grad2z(interp2z(w) * interp2z(w)) + visc * (grad2z(grad2z(w)))
            + alpha*interp2z(s - s_ref) )
    end
end
