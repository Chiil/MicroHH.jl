## Load the packages.
using LoopVectorization
using Interpolations
using Statistics
using BenchmarkTools
using PyPlot


## Create the refine function.
function refine_field_int!(hi, hi_tmp, lo, n_hi, n_lo)
    if n_hi ÷ n_lo == 2
        @tturbo for k in 0:size(lo, 3)-3
            for j in 0:size(lo, 2)-3
                for i in 0:size(lo, 1)-3
                    i_lo = i + 2
                    j_lo = j + 2
                    k_lo = k + 2
                    i_hi = i*2 + 2
                    j_hi = j*2 + 2
                    k_hi = k*2 + 2
                    for kk in 0:1, jj in 0:1, ii in 0:1
                        #hi_tmp[i_hi+ii, j_hi+jj, k_hi+kk] = lo[i_lo, j_lo, k_lo]
                        hi[i_hi+ii, j_hi+jj, k_hi+kk] = lo[i_lo, j_lo, k_lo]
                    end
                end
            end
        end

        #hi[:, :, :] = hi_tmp[:, :, :]

        # hi_tmp[1, :] = hi_tmp[end-1, :]
        # hi_tmp[end, :] = hi_tmp[2, :]
        # hi_tmp[:, 1] = hi_tmp[:, end-1]
        # hi_tmp[:, end] = hi_tmp[:, 2]

        # @tturbo for j in 2:size(hi, 2)-1
        #     for i in 2:size(hi, 1)-1
        #         hi[i, j] = ( 1//16*hi_tmp[i-1, j-1] + 1//8*hi_tmp[i, j-1] + 1//16*hi_tmp[i+1, j-1]
        #                    + 1// 8*hi_tmp[i-1, j  ] + 1//4*hi_tmp[i, j  ] + 1// 8*hi_tmp[i+1, j  ]
        #                    + 1//16*hi_tmp[i-1, j+1] + 1//8*hi_tmp[i, j+1] + 1//16*hi_tmp[i+1, j+1] )
        #     end
        # end
    elseif n_hi ÷ n_lo == 3
        @tturbo for k in 0:size(lo, 3)-3
            for j in 0:size(lo, 2)-3
                for i in 0:size(lo, 1)-3
                    i_lo = i + 2
                    j_lo = j + 2
                    k_lo = k + 2
                    i_hi = i*3 + 2
                    j_hi = j*3 + 2
                    k_hi = k*3 + 2
                    for kk in 0:2, jj in 0:2, ii in 0:2
                        #hi_tmp[i_hi+ii, j_hi+jj, k_hi+kk] = lo[i_lo, j_lo, k_lo]
                        hi[i_hi+ii, j_hi+jj, k_hi+kk] = lo[i_lo, j_lo, k_lo]
                    end
                end
            end
        end

        #hi[:, :, :] = hi_tmp[:, :, :]

        # hi_tmp[1, :] = hi_tmp[end-1, :]
        # hi_tmp[end, :] = hi_tmp[2, :]
        # hi_tmp[:, 1] = hi_tmp[:, end-1]
        # hi_tmp[:, end] = hi_tmp[:, 2]

        # @tturbo for j in 2:size(hi, 2)-1
        #     for i in 2:size(hi, 1)-1
        #         hi[i, j] = ( 1//9*hi_tmp[i-1, j-1] + 1//9*hi_tmp[i, j-1] + 1//9*hi_tmp[i+1, j-1]
        #                    + 1//9*hi_tmp[i-1, j  ] + 1//9*hi_tmp[i, j  ] + 1//9*hi_tmp[i+1, j  ]
        #                    + 1//9*hi_tmp[i-1, j+1] + 1//9*hi_tmp[i, j+1] + 1//9*hi_tmp[i+1, j+1] )
        #     end
        # end
    else
        throw(DomainError(n_hi/n_lo, "Refinement should be 2 or 3"))
    end
end


function refine_field_int_ref!(hi, hi_tmp, lo, x_hi, x_lo, n_hi, n_lo)
    interp = LinearInterpolation((x_lo, x_lo, x_lo), lo)
    hi[2:end-1, 2:end-1, 2:end-1] = interp(x_hi[2:end-1], x_hi[2:end-1], x_hi[2:end-1])
end

## Set up the grids.
n_hi = 384
n_lo = n_hi ÷ 3

a_lo = rand(n_lo, n_lo, n_lo)
a_hi = zeros(n_hi, n_hi, n_hi)

a_hi_gc = zeros(n_hi+2, n_hi+2, n_hi+2)
a_lo_gc = zeros(n_lo+2, n_lo+2, n_lo+2)

a_hi_gc[2:n_hi+1, 2:n_hi+1, 2:n_hi+1] = a_hi[:, :, :]
a_lo_gc[2:n_lo+1, 2:n_lo+1, 2:n_lo+1] = a_lo[:, :, :]

# Cyclic BCs in all 3 directions.
a_hi_gc[1, :, :] = a_hi_gc[end-1, :, :]
a_hi_gc[end, :, :] = a_hi_gc[2, :, :]
a_hi_gc[:, 1, :] = a_hi_gc[:, end-1, :]
a_hi_gc[:, end, :] = a_hi_gc[:, 2, :]
a_hi_gc[:, :, 1] = a_hi_gc[:, :, end-1, :]
a_hi_gc[:, :, end] = a_hi_gc[:, :, 2]

a_lo_gc[1, :, :] = a_lo_gc[end-1, :, :]
a_lo_gc[end, :, :] = a_lo_gc[2, :, :]
a_lo_gc[:, 1, :] = a_lo_gc[:, end-1, :]
a_lo_gc[:, end, :] = a_lo_gc[:, 2, :]
a_lo_gc[:, :, 1] = a_lo_gc[:, :, end-1, :]
a_lo_gc[:, :, end] = a_lo_gc[:, :, 2]


dx_hi = 1/n_hi; dx_lo = 1/n_lo
x_hi = dx_hi/2:dx_hi:1 |> collect
x_lo = dx_lo/2:dx_lo:1 |> collect
xh_hi = 0:dx_hi:1 |> collect
xh_lo = 0:dx_lo:1 |> collect

x_hi_gc = -dx_hi/2:dx_hi:1+dx_hi |> collect
x_lo_gc = -dx_lo/2:dx_lo:1+dx_lo |> collect

## Compute.
a_hi_tmp = zeros(size(a_hi_gc))
@btime refine_field_int!(a_hi_gc, a_hi_tmp, a_lo_gc, n_hi, n_lo)
a_hi_int = a_hi_gc[2:end-1, 2:end-1, 2:end-1]

## Compute a reference using Interpolations.jl
a_hi_tmp = zeros(size(a_hi_gc))
@btime refine_field_int_ref!(a_hi_gc, a_hi_tmp, a_lo_gc, x_hi_gc, x_lo_gc, n_hi, n_lo)
a_hi_int_ref = a_hi_gc[2:end-1, 2:end-1, 2:end-1]

println("Mean equal to lo: ", mean(a_lo) ≈ mean(a_hi_int))
println("Mean equal to ref: ", mean(a_hi_int_ref) ≈ mean(a_hi_int))
println("Interpolations equal: ", a_hi_int ≈ a_hi_int_ref)


## Plot the output.
figure(figsize=(10, 5))
subplot(121)
pcolormesh(xh_lo, xh_lo, a_lo[:, :, 1], vmin=0, vmax=1)
subplot(122)
pcolormesh(xh_hi, xh_hi, a_hi_int[:, :, 1], vmin=0, vmax=1)
tight_layout()
display(gcf())
show()