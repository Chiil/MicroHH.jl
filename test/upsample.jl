## Load the packages.
using LoopVectorization
using Interpolations
using OffsetArrays
using Statistics
using BenchmarkTools
using PyPlot
using Printf


## Define the upsample functions.
function upsample_nn!(hi, lo, itot, jtot, ktot, ifac, jfac, kfac, ioff, joff, koff)
    @tturbo for k in 0:ktot-1, j in 0:jtot-1, i in 0:itot-1
        i_lo = i + 2; j_lo = j + 2; k_lo = k + 2
        i_hi = ifac*i + 2; j_hi = jfac*j + 2; k_hi = kfac*k + 2
        for kk in 0:kfac-1, jj in 0:jfac-1, ii in 0:ifac-1
            iii = floor(Int, ioff + 1/ifac * ii)
            jjj = floor(Int, joff + 1/jfac * jj)
            kkk = floor(Int, koff + 1/kfac * kk)
            hi[i_hi+ii, j_hi+jj, k_hi+kk] = lo[i_lo+iii, j_lo+jjj, k_lo+kkk]
        end
    end
end

function upsample_nn_ref!(hi, lo, x_hi, y_hi, z_hi, x_lo, y_lo, z_lo)
    interp = interpolate((x_lo, y_lo, z_lo), lo, (Gridded(Constant()), Gridded(Constant()), Gridded(Constant())))
    hi[2:end-1, 2:end-1, 2:end-1] .= interp(x_hi[2:end-1], y_hi[2:end-1], z_hi[2:end-1])
end

macro upsample_lin_fast(suffix, ifac, jfac, kfac, ioff, joff, koff)
    ex_inner = []
    for kk in 0:kfac-1, jj in 0:jfac-1, ii in 0:ifac-1
        push!(ex_inner, :(i_hi = $ifac*i+$ii+is_hi), :(j_hi = $jfac*j+$jj+js_hi), :(k_hi = $kfac*k+$kk+ks_hi))

        i_pos = -1/2 + ioff + (1/2 - ioff + ii)/ifac
        j_pos = -1/2 + joff + (1/2 - joff + jj)/jfac
        k_pos = -1/2 + koff + (1/2 - koff + kk)/kfac

        i_lo_off = floor(Int, i_pos) + is_lo
        j_lo_off = floor(Int, j_pos) + js_lo
        k_lo_off = floor(Int, k_pos) + ks_lo
        push!(ex_inner, :(i_lo = i + $i_lo_off), :(j_lo = j + $j_lo_off), :(k_lo = k + $k_lo_off))

        fi = mod(i_pos, 1)
        fj = mod(j_pos, 1)
        fk = mod(k_pos, 1)

        coef = zeros(2, 2, 2)
        coef[1, 1, 1] = (((1-fi)*fj+fi-1)*fk+(fi-1)*fj-fi+1)
        coef[2, 1, 1] = ((fi*fj-fi)*fk-fi*fj+fi)
        coef[1, 2, 1] = ((fi-1)*fj*fk+(1-fi)*fj)
        coef[2, 2, 1] = (fi*fj-fi*fj*fk)
        coef[1, 1, 2] = ((fi-1)*fj-fi+1)*fk
        coef[2, 1, 2] = (fi-fi*fj)*fk
        coef[1, 2, 2] = (1-fi)*fj*fk
        coef[2, 2, 2] = fi*fj*fk

        rhs_list = []
        for kc in 1:2, jc in 1:2, ic in 1:2
            if !(coef[ic, jc, kc] ≈ 0)
                push!(rhs_list, :( $(coef[ic, jc, kc])*lo[i_lo+$(ic-1), j_lo+$(jc-1), k_lo+$(kc-1)] ))
            end
        end

        rhs = Expr(:call, :+, rhs_list...)

        push!(ex_inner, :(@inbounds hi[i_hi, j_hi, k_hi] = $rhs))
    end

    ex_inner_block = Expr(:block, ex_inner...)

    name = Symbol(@sprintf "upsample_lin_%s!" suffix)
    ex = quote
        function $name(hi, lo, itot, jtot, ktot, is_hi, js_hi, ks_hi, is_lo, js_lo, ks_lo)
            Threads.@threads for k in 0:ktot-1
                for j in 0:jtot-1
                    for i in 0:itot-1
                        $ex_inner_block
                    end
                end
            end
        end
    end
    return esc(ex)
end

function upsample_lin!(
        hi, lo,
        itot, jtot, ktot,
        is_hi, js_hi, ks_hi, is_lo, js_lo, ks_lo,
        ifac, jfac, kfac, ioff, joff, koff)

    Threads.@threads for k in 0:ktot-1
        for j in 0:jtot-1
            for i in 0:itot-1
                for kk in 0:kfac-1, jj in 0:jfac-1, ii in 0:ifac-1
                    i_hi = ifac*i+ii+2; j_hi = jfac*j+jj+2; k_hi = kfac*k+kk+2;

                    # Determine fractional distance and west, south, and bot point.
                    i_pos = -1/2 + ioff + (1/2 - ioff + ii)/ifac
                    j_pos = -1/2 + joff + (1/2 - joff + jj)/jfac
                    k_pos = -1/2 + koff + (1/2 - koff + kk)/kfac

                    fi = mod(i_pos, 1)
                    fj = mod(j_pos, 1)
                    fk = mod(k_pos, 1)

                    coef_111 = (((1-fi)*fj+fi-1)*fk+(fi-1)*fj-fi+1)
                    coef_211 = ((fi*fj-fi)*fk-fi*fj+fi)
                    coef_121 = ((fi-1)*fj*fk+(1-fi)*fj)
                    coef_221 = (fi*fj-fi*fj*fk)
                    coef_112 = ((fi-1)*fj-fi+1)*fk
                    coef_212 = (fi-fi*fj)*fk
                    coef_122 = (1-fi)*fj*fk
                    coef_222 = fi*fj*fk

                    i_lo = i + floor(Int, i_pos) + 2
                    j_lo = j + floor(Int, j_pos) + 2
                    k_lo = k + floor(Int, k_pos) + 2

                    @inbounds hi[i_hi, j_hi, k_hi] = (
                        + coef_111 * lo[i_lo  , j_lo  , k_lo  ]
                        + coef_211 * lo[i_lo+1, j_lo  , k_lo  ]
                        + coef_121 * lo[i_lo  , j_lo+1, k_lo  ]
                        + coef_221 * lo[i_lo+1, j_lo+1, k_lo  ]
                        + coef_112 * lo[i_lo  , j_lo  , k_lo+1]
                        + coef_212 * lo[i_lo+1, j_lo  , k_lo+1]
                        + coef_122 * lo[i_lo  , j_lo+1, k_lo+1]
                        + coef_222 * lo[i_lo+1, j_lo+1, k_lo+1] )
                end
            end
        end
    end
end

function upsample_lin_ref!(hi, lo, x_hi, y_hi, z_hi, x_lo, y_lo, z_lo)
    interp = interpolate((x_lo, y_lo, z_lo), lo, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    hi[2:end-1, 2:end-1, 2:end-1] .= interp(x_hi[2:end-1], y_hi[2:end-1], z_hi[2:end-1])
end


## Set up the grids.
itot_lo = 256; jtot_lo = 192; ktot_lo = 128
is_lo = 2; js_lo = 2; ks_lo = 2
is_hi = 2; js_hi = 2; ks_hi = 2
ifac = 2; jfac = 2; kfac = 2

a_lo = rand(itot_lo + 2, jtot_lo + 2, ktot_lo + 2)
a_hi_ref = zeros(itot_lo*ifac + 2, jtot_lo*jfac + 2, ktot_lo*kfac + 2)
a_hi = zeros(itot_lo*ifac + 2, jtot_lo*jfac + 2, ktot_lo*kfac + 2)
a_hi_fast = zeros(itot_lo*ifac + 2, jtot_lo*jfac + 2, ktot_lo*kfac + 2)

dx_lo = 1/itot_lo; dy_lo = 1/jtot_lo; dz_lo = 1/ktot_lo
x_lo = -dx_lo/2:dx_lo:1+dx_lo |> collect
y_lo = -dy_lo/2:dy_lo:1+dy_lo |> collect
z_lo = -dz_lo/2:dz_lo:1+dz_lo |> collect
xh_lo = -dx_lo:dx_lo:1 |> collect
yh_lo = -dy_lo:dy_lo:1 |> collect

dx_hi = dx_lo/ifac; dy_hi = dy_lo/jfac; dz_hi = dz_lo/kfac
x_hi = -dx_hi/2:dx_hi:1+dx_hi |> collect
y_hi = -dy_hi/2:dy_hi:1+dy_hi |> collect
z_hi = -dz_hi/2:dz_hi:1+dz_hi |> collect
xh_hi = -dx_hi:dx_hi:1 |> collect
yh_hi = -dy_hi:dy_hi:1 |> collect


## Compute a reference using Interpolations.jl
@btime upsample_lin_ref!(a_hi_ref, a_lo, x_hi, y_hi, z_hi, x_lo, y_lo, z_lo)
@btime upsample_lin!(a_hi, a_lo, itot_lo, jtot_lo, ktot_lo,
    is_hi, js_hi, ks_hi, is_lo, js_lo, ks_lo,
    ifac, jfac, kfac, 0, 0, 0)
@upsample_lin_fast("222", 2, 2, 2, 0, 0, 0)
@btime upsample_lin_222!(
    a_hi_fast, a_lo, itot_lo, jtot_lo, ktot_lo,
    is_hi, js_hi, ks_hi, is_lo, js_lo, ks_lo)

a_lo_int = @view a_lo[2:end-1, 2:end-1, 2:end-1]
a_hi_ref_int = @view a_hi_ref[2:end-1, 2:end-1, 2:end-1]
a_hi_int = @view a_hi[2:end-1, 2:end-1, 2:end-1]
a_hi_fast_int = @view a_hi_fast[2:end-1, 2:end-1, 2:end-1]

println("Standard equal to ref: ", a_hi_int ≈ a_hi_ref_int)
println("Fast equal to ref: ", a_hi_fast_int ≈ a_hi_ref_int)


## Plot the output.
x_lo_int = @view x_lo[2:end-1]
x_hi_int = @view x_hi[2:end-1]

figure()
plot(x_hi_int, a_hi_int[:, 1, 1], "C0-o")
plot(x_hi_int, a_hi_ref_int[:, 1, 1], "C1-^")
plot(x_hi_int, a_hi_fast_int[:, 1, 1], "C2-*")
plot(x_lo_int, a_lo_int[:, 1, 1], "k:+")
tight_layout()
display(gcf())


## Compute a reference using Interpolations.jl
u_lo = rand(itot_lo + 2, jtot_lo + 2, ktot_lo + 2)
u_hi_ref = zeros(itot_lo*ifac + 2, jtot_lo*jfac + 2, ktot_lo*kfac + 2)
u_hi = zeros(itot_lo*ifac + 2, jtot_lo*jfac + 2, ktot_lo*kfac + 2)
u_hi_fast = zeros(itot_lo*ifac + 2, jtot_lo*jfac + 2, ktot_lo*kfac + 2)

@btime upsample_lin_ref!(u_hi_ref, u_lo, xh_hi, y_hi, z_hi, xh_lo, y_lo, z_lo)
@btime upsample_lin!(
    u_hi, u_lo, itot_lo, jtot_lo, ktot_lo,
    is_hi, js_hi, ks_hi, is_lo, js_lo, ks_lo,
    ifac, jfac, kfac, 0.5, 0, 0)
@upsample_lin_fast("222_u", 2, 2, 2, 0.5, 0, 0)
@btime upsample_lin_222_u!(
    u_hi_fast, u_lo, itot_lo, jtot_lo, ktot_lo,
    is_hi, js_hi, ks_hi, is_lo, js_lo, ks_lo)

u_lo_int = @view u_lo[2:end-1, 2:end-1, 2:end-1]
u_hi_ref_int = @view u_hi_ref[2:end-1, 2:end-1, 2:end-1]
u_hi_int = @view u_hi[2:end-1, 2:end-1, 2:end-1]
u_hi_fast_int = @view u_hi_fast[2:end-1, 2:end-1, 2:end-1]

println("Standard equal to ref: ", u_hi_int ≈ u_hi_ref_int)
println("Fast equal to ref: ", u_hi_fast_int ≈ u_hi_ref_int)


## Plot the output.
xh_lo_int = @view xh_lo[2:end-1]
xh_hi_int = @view xh_hi[2:end-1]

figure()
plot(xh_hi_int, u_hi_int[:, 1, 1], "C0-o")
plot(xh_hi_int, u_hi_ref_int[:, 1, 1], "C1-^")
plot(xh_hi_int, u_hi_fast_int[:, 1, 1], "C2-*")
plot(xh_lo_int, u_lo_int[:, 1, 1], "k:+")
tight_layout()
display(gcf())

show()
