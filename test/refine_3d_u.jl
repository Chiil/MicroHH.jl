## Load the packages.
using LoopVectorization
using Interpolations
using OffsetArrays
using Statistics
using BenchmarkTools
using PyPlot


## Create the refine functions.
function refine_field_int!(hi, hi_tmp, lo, n_hi, n_lo)
    if n_hi ÷ n_lo == 2
        @tturbo for k in 0:size(lo, 3)-3, j in 0:size(lo, 2)-3, i in 0:size(lo, 1)-3
            i_lo = i + 2
            j_lo = j + 2
            k_lo = k + 2
            i_hi = i*2 + 2
            j_hi = j*2 + 2
            k_hi = k*2 + 2
            for kk in 0:1, jj in 0:1
                hi_tmp[i_hi, j_hi+jj, k_hi+kk] = lo[i_lo, j_lo, k_lo]
                hi_tmp[i_hi+1, j_hi+jj, k_hi+kk] = 1//2*(lo[i_lo, j_lo, k_lo] + lo[i_lo+1, j_lo, k_lo])
            end
        end

        @tturbo a_hi_tmp[1, :, :] = a_hi_tmp[end-1, :, :]
        @tturbo a_hi_tmp[end, :, :] = a_hi_tmp[2, :, :]
        @tturbo a_hi_tmp[:, 1, :] = a_hi_tmp[:, end-1, :]
        @tturbo a_hi_tmp[:, end, :] = a_hi_tmp[:, 2, :]
        @tturbo a_hi_tmp[:, :, 1] = a_hi_tmp[:, :, end-1]
        @tturbo a_hi_tmp[:, :, end] = a_hi_tmp[:, :, 2]

        coef = OffsetArray(zeros(3, 3), -1:1, -1:1)
        for kk in -1:1, jj in -1:1
            coef[jj, kk] = (1//2)^(2 + abs(jj) + abs(kk))
        end

        @tturbo for k in 2:size(hi, 3)-1, j in 2:size(hi, 2)-1, i in 2:size(hi, 1)-1
            hi[i, j, k] = 0
            for kk in -1:1, jj in -1:1
                hi[i, j, k] += coef[jj, kk] * hi_tmp[i, j+jj, k+kk]
            end
        end
    elseif n_hi ÷ n_lo == 3
        @tturbo for k in 0:size(lo, 3)-3, j in 0:size(lo, 2)-3, i in 0:size(lo, 1)-3
            i_lo = i + 2
            j_lo = j + 2
            k_lo = k + 2
            i_hi = i*3 + 2
            j_hi = j*3 + 2
            k_hi = k*3 + 2
            for kk in 0:2, jj in 0:2
                hi_tmp[i_hi  , j_hi+jj, k_hi+kk] = lo[i_lo  , j_lo, k_lo]
                hi_tmp[i_hi+1, j_hi+jj, k_hi+kk] = lo[i_lo  , j_lo, k_lo]
                hi_tmp[i_hi+2, j_hi+jj, k_hi+kk] = lo[i_lo+1, j_lo, k_lo]
            end
        end

        @tturbo a_hi_tmp[1, :, :] = a_hi_tmp[end-1, :, :]
        @tturbo a_hi_tmp[end, :, :] = a_hi_tmp[2, :, :]
        @tturbo a_hi_tmp[:, 1, :] = a_hi_tmp[:, end-1, :]
        @tturbo a_hi_tmp[:, end, :] = a_hi_tmp[:, 2, :]
        @tturbo a_hi_tmp[:, :, 1] = a_hi_tmp[:, :, end-1]
        @tturbo a_hi_tmp[:, :, end] = a_hi_tmp[:, :, 2]

        @tturbo for k in 2:size(hi, 3)-1, j in 2:size(hi, 2)-1, i in 2:size(hi, 1)-1
            hi[i, j, k] = 0
            for kk in -1:1, jj in -1:1, ii in -1:1
                hi[i, j, k] += 1//27*hi_tmp[i+ii, j+jj, k+kk]
            end
        end
    else
        throw(DomainError(n_hi/n_lo, "Refinement should be 2 or 3"))
    end
end

function refine_field_int_ref!(hi, hi_tmp, lo, x_hi, xh_hi, x_lo, xh_lo, n_hi, n_lo)
    interp = interpolate((xh_lo, x_lo, x_lo), lo, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    # interp = interpolate((xh_lo, x_lo, x_lo), lo, (Gridded(Constant()), Gridded(Constant()), Gridded(Constant())))
    hi[2:end-1, 2:end-1, 2:end-1] = interp(xh_hi[2:end-1], x_hi[2:end-1], x_hi[2:end-1])
end


## Set up the grids.
n_hi = 384
n_lo = n_hi ÷ 2

a_lo = rand(n_lo, n_lo, n_lo)
a_hi = zeros(n_hi, n_hi, n_hi)

a_hi_gc = zeros(n_hi+2, n_hi+2, n_hi+2)
a_lo_gc = zeros(n_lo+2, n_lo+2, n_lo+2)

a_hi_gc[2:n_hi+1, 2:n_hi+1, 2:n_hi+1] = a_hi[:, :, :]
a_lo_gc[2:n_lo+1, 2:n_lo+1, 2:n_lo+1] = a_lo[:, :, :]


## Cyclic BCs in all 3 directions.
a_hi_gc[1, :, :] = a_hi_gc[end-1, :, :]
a_hi_gc[end, :, :] = a_hi_gc[2, :, :]
a_hi_gc[:, 1, :] = a_hi_gc[:, end-1, :]
a_hi_gc[:, end, :] = a_hi_gc[:, 2, :]
a_hi_gc[:, :, 1] = a_hi_gc[:, :, end-1]
a_hi_gc[:, :, end] = a_hi_gc[:, :, 2]

a_lo_gc[1, :, :] = a_lo_gc[end-1, :, :]
a_lo_gc[end, :, :] = a_lo_gc[2, :, :]
a_lo_gc[:, 1, :] = a_lo_gc[:, end-1, :]
a_lo_gc[:, end, :] = a_lo_gc[:, 2, :]
a_lo_gc[:, :, 1] = a_lo_gc[:, :, end-1]
a_lo_gc[:, :, end] = a_lo_gc[:, :, 2]


dx_hi = 1/n_hi; dx_lo = 1/n_lo
x_hi = dx_hi/2:dx_hi:1 |> collect
x_lo = dx_lo/2:dx_lo:1 |> collect
xh_hi = 0:dx_hi:1-dx_hi/2 |> collect
xh_lo = 0:dx_lo:1-dx_lo/2 |> collect

x_hi_gc = -dx_hi/2:dx_hi:1+dx_hi |> collect
x_lo_gc = -dx_lo/2:dx_lo:1+dx_lo |> collect
xh_hi_gc = -dx_hi:dx_hi:1+dx_hi/2 |> collect
xh_lo_gc = -dx_lo:dx_lo:1+dx_lo/2 |> collect


## Compute.
a_hi_tmp = zeros(size(a_hi_gc))
@btime refine_field_int!(a_hi_gc, a_hi_tmp, a_lo_gc, n_hi, n_lo)
a_hi_int = a_hi_gc[2:end-1, 2:end-1, 2:end-1]


## Compute a reference using Interpolations.jl
a_hi_tmp = zeros(size(a_hi_gc))
a_hi_gc[:, :, :] .= 0.
@btime refine_field_int_ref!(a_hi_gc, a_hi_tmp, a_lo_gc, x_hi_gc, xh_hi_gc, x_lo_gc, xh_lo_gc, n_hi, n_lo)
a_hi_int_ref = a_hi_gc[2:end-1, 2:end-1, 2:end-1]

println("Mean equal to lo: ", mean(a_lo) ≈ mean(a_hi_int))
println("Mean equal to ref: ", mean(a_hi_int_ref) ≈ mean(a_hi_int))
println("Interpolations equal: ", a_hi_int ≈ a_hi_int_ref)

x_hi_p = copy(x_hi)
x_lo_p = copy(x_lo)
pushfirst!(x_hi_p, 0)
pushfirst!(x_lo_p, 0)
xh_hi_p = copy(xh_hi)
xh_lo_p = copy(xh_lo)
push!(xh_hi_p, 1)
push!(xh_lo_p, 1)


## Plot the output.
figure(figsize=(10, 5))
subplot(121)
pcolormesh(x_lo_p, xh_lo_p, a_lo[:, :, 1]', vmin=0, vmax=1)
xlim(0, 1)
ylim(0, 1)
subplot(122)
pcolormesh(x_hi_p, xh_hi_p, a_hi_int[:, :, 1]', vmin=0, vmax=1)
xlim(0, 1)
ylim(0, 1)
tight_layout()
display(gcf())
show()
