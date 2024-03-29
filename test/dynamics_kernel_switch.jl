## Packages and settings.
include("../src/StencilBuilder.jl")

using BenchmarkTools
using LoopVectorization
using .StencilBuilder

float_type = Float32


## CPU dynamics kernel.
macro make_kernel(do_advec::Bool, do_diff::Bool, do_buoy::Bool)
    ex_rhs_l = []
    if do_advec
        push!(ex_rhs_l, :( -grad2x(interp2z(u) * interp2x(w)) ))
        push!(ex_rhs_l, :( -grad2y(interp2z(v) * interp2y(w)) ))
        push!(ex_rhs_l, :( -grad2z(interp2z(w) * interp2z(w)) ))
    end
    if do_diff
        push!(ex_rhs_l, :( visc * (grad2x(grad2x(w))) ))
        push!(ex_rhs_l, :( visc * (grad2y(grad2y(w))) ))
        push!(ex_rhs_l, :( visc * (grad2z(grad2z(w))) ))
    end
    if do_buoy
        push!(ex_rhs_l, :( alpha*interp2z(s) ))
    end

    if length(ex_rhs_l) > 0
        ex_rhs = Expr(:call, :+, ex_rhs_l...)

        ex = quote
            function dynamics_w_kernel!(
                    wt, u, v, w, s,
                    visc, alpha,
                    dxi, dyi, dzi, dzhi,
                    is, ie, js, je, ks, ke)

                @fast3d begin
                    @fd (wt, u, v, w, s) wt += $ex_rhs
                end
            end
        end
    else
        ex = quote
            function dynamics_w_kernel!(
                    wt, u, v, w, s,
                    visc, alpha,
                    dxi, dyi, dzi, dzhi,
                    is, ie, js, je, ks, ke)
            end
        end
    end

    return esc(ex)
end


## User input
itot = 384; jtot = 384; ktot = 384;
igc = 1; jgc = 1; kgc = 1;

do_advec = true; do_diff = true; do_buoy = true

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
visc = convert(float_type, 1.)
alpha = convert(float_type, 9.81/300)


## Run CPU kernel
@eval @make_kernel($do_advec, $do_diff, $do_buoy)

@btime dynamics_w_kernel!(
    $wt, $u, $v, $w, $s,
    $visc, $alpha,
    $dxi, $dyi, $dzi, $dzhi,
    $is, $ie, $js, $je, $ks, $ke)

