## Load the packages.
using LoopVectorization
using Interpolations
using Statistics
using BenchmarkTools
using PyPlot
using Printf


## Define the upsample functions.
function make_inner_loop!(ex_list, ii, jj, kk, ifac, jfac, kfac, ioff, joff, koff, is_top)
    push!(ex_list, :(i_hi = $ifac*i+$ii+is_hi), :(j_hi = $jfac*j+$jj+js_hi))
    if is_top
        push!(ex_list, :(k_hi = ke_hi))
    else
        push!(ex_list, :(k_hi = $kfac*k+$kk+ks_hi))
    end

    i_pos = -1/2 + ioff + (1/2 - ioff + ii)/ifac
    j_pos = -1/2 + joff + (1/2 - joff + jj)/jfac
    k_pos = -1/2 + koff + (1/2 - koff + kk)/kfac

    i_lo_off = floor(Int, i_pos) + is_lo
    j_lo_off = floor(Int, j_pos) + js_lo
    k_lo_off = floor(Int, k_pos) + ks_lo

    if is_top
        push!(ex_list, :(i_lo = i + $i_lo_off), :(j_lo = j + $j_lo_off), :(k_lo = ke_lo))
    else
        push!(ex_list, :(i_lo = i + $i_lo_off), :(j_lo = j + $j_lo_off), :(k_lo = k + $k_lo_off))
    end

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

    kc_end = is_top ? 1 : 2

    rhs_list = []
    for kc in 1:kc_end, jc in 1:2, ic in 1:2
        if !(coef[ic, jc, kc] ≈ 0)
            push!(rhs_list, :( $(coef[ic, jc, kc])*lo[i_lo+$(ic-1), j_lo+$(jc-1), k_lo+$(kc-1)] ))
        end
    end

    rhs = Expr(:call, :+, rhs_list...)

    push!(ex_list, :(hi[i_hi, j_hi, k_hi] = $rhs))
end

macro upsample_lin_fast(suffix, ifac, jfac, kfac, ioff, joff, koff, add_top)
    ex_inner = []
    for kk in 0:kfac-1, jj in 0:jfac-1, ii in 0:ifac-1
        make_inner_loop!(ex_inner, ii, jj, kk, ifac, jfac, kfac, ioff, joff, koff, false)
    end
    ex_inner_block = quote
        Threads.@threads for k in 0:ktot-1
            for j in 0:jtot-1
                @inbounds @simd for i in 0:itot-1
                    $(Expr(:block, ex_inner...))
                end
            end
        end
    end

    if add_top
        ex_top = []
        kk = 0
        for jj in 0:jfac-1, ii in 0:ifac-1
            make_inner_loop!(ex_top, ii, jj, kk, ifac, jfac, kfac, ioff, joff, koff, true)
        end
        ex_top_block = quote
            Threads.@threads for j in 0:jtot-1
                @inbounds @simd for i in 0:itot-1
                    $(Expr(:block, ex_top...))
                end
            end
        end
    else
        ex_top_block = :()
    end

    name = Symbol(@sprintf "upsample_lin_%s!" suffix)
    ex = quote
        function $name(hi, lo, itot, jtot, ktot, is_hi, js_hi, ks_hi, is_lo, js_lo, ks_lo)
            $ex_inner_block
            $ex_top_block
        end
    end
    println(ex)
    return esc(ex)
end

function upsample_lin_ref!(hi, lo, x_hi, y_hi, z_hi, x_lo, y_lo, z_lo, is_hi, js_hi, ks_hi, ie_hi, je_hi, ke_hi)
    interp = interpolate((x_lo, y_lo, z_lo), lo, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    hi[is_hi:ie_hi, js_hi:je_hi, ks_hi:ke_hi] .= interp(x_hi[is_hi:ie_hi], y_hi[js_hi:je_hi], z_hi[ks_hi:ke_hi])
end


## Set up the grids.
itot_lo = 256; jtot_lo = 192; ktot_lo = 128
itot_hi = 512; jtot_hi = 384; ktot_hi = 256
igc_lo = 1; jgc_lo = 1; kgc_lo = 1
igc_hi = 1; jgc_hi = 1; kgc_hi = 1

# Derive dimensions
ifac = itot_hi÷itot_lo; jfac = jtot_hi÷jtot_lo; kfac = ktot_hi÷ktot_lo;
icells_lo = itot_lo + 2igc_lo; jcells_lo = jtot_lo + 2jgc_lo; kcells_lo = ktot_lo + 2kgc_lo
icells_hi = itot_hi + 2igc_hi; jcells_hi = jtot_hi + 2jgc_hi; kcells_hi = ktot_hi + 2kgc_hi
is_lo = igc_lo+1; js_lo = jgc_lo+1; ks_lo = kgc_lo+1
is_hi = igc_hi+1; js_hi = jgc_hi+1; ks_hi = kgc_hi+1
ie_lo = igc_lo+itot_lo; je_lo = jgc_lo+jtot_lo; ke_lo = kgc_lo+ktot_lo+1
ie_hi = igc_hi+itot_hi; je_hi = jgc_hi+jtot_hi; ke_hi = kgc_hi+ktot_hi+1

w_lo = rand(icells_lo, jcells_lo, kcells_lo)
w_hi_ref = zeros(icells_hi, jcells_hi, kcells_hi)
w_hi_fast = zeros(icells_hi, jcells_hi, kcells_hi)

dx_lo = 1/itot_lo; dy_lo = 1/jtot_lo; dz_lo = 1/ktot_lo
x_lo = collect(range(0.5-igc_lo, length=icells_lo)) .* dx_lo
y_lo = collect(range(0.5-jgc_lo, length=jcells_lo)) .* dy_lo
z_lo = collect(range(0.5-kgc_lo, length=kcells_lo)) .* dz_lo
xh_lo = collect(range(-igc_lo, length=icells_lo)) .* dx_lo
yh_lo = collect(range(-jgc_lo, length=jcells_lo)) .* dy_lo
zh_lo = collect(range(-kgc_lo, length=kcells_lo)) .* dz_lo

dx_hi = 1/itot_hi; dy_hi = 1/jtot_hi; dz_hi = 1/ktot_hi
x_hi = collect(range(0.5-igc_hi, length=icells_hi)) .* dx_hi
y_hi = collect(range(0.5-jgc_hi, length=jcells_hi)) .* dy_hi
z_hi = collect(range(0.5-kgc_hi, length=kcells_hi)) .* dz_hi
xh_hi = collect(range(-igc_hi, length=icells_hi)) .* dx_hi
yh_hi = collect(range(-jgc_hi, length=jcells_hi)) .* dy_hi
zh_hi = collect(range(-kgc_hi, length=kcells_hi)) .* dz_hi


## Compute a reference using Interpolations.jl
@upsample_lin_fast("222", 2, 2, 2, 0, 0, 0.5, true)
@btime upsample_lin_ref!(
    w_hi_ref, w_lo, x_hi, y_hi, zh_hi, x_lo, y_lo, zh_lo,
    is_hi, js_hi, ks_hi, ie_hi, je_hi, ke_hi)
@btime upsample_lin_222!(
    w_hi_fast, w_lo, itot_lo, jtot_lo, ktot_lo,
    is_hi, js_hi, ks_hi, is_lo, js_lo, ks_lo)

w_lo_int = @view w_lo[is_lo:ie_lo, js_lo:je_lo, ks_lo:ke_lo]
w_hi_ref_int = @view w_hi_ref[is_hi:ie_hi, js_hi:je_hi, ks_hi:ke_hi]
w_hi_fast_int = @view w_hi_fast[is_hi:ie_hi, js_hi:je_hi, ks_hi:ke_hi]

println("Fast equal to ref: ", w_hi_fast_int ≈ w_hi_ref_int)


## Plot the output.
zh_lo_int = @view zh_lo[ks_lo:ke_lo]
zh_hi_int = @view zh_hi[ks_hi:ke_hi]

figure()
plot(zh_hi_int, w_hi_ref_int[1, 1, :], "C1-^")
plot(zh_hi_int, w_hi_fast_int[1, 1, :], "C2-*")
plot(zh_lo_int, w_lo_int[1, 1, :], "k:+")
tight_layout()
display(gcf())
show()
