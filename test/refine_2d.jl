## Load the packages.
using LoopVectorization
using Statistics
using BenchmarkTools
using PyPlot


## Create the refine function.
function refine_field_int!(hi, hi_tmp, lo, n_hi, n_lo)
    if n_hi ÷ n_lo == 2
        @tturbo for j in 0:size(lo, 2)-3
            for i in 0:size(lo, 1)-3
                i_lo = i + 2
                j_lo = j + 2
                i_hi = i*2 + 2
                j_hi = j*2 + 2
                hi_tmp[i_hi  , j_hi  ] = lo[i_lo, j_lo]
                hi_tmp[i_hi+1, j_hi  ] = lo[i_lo, j_lo]
                hi_tmp[i_hi  , j_hi+1] = lo[i_lo, j_lo]
                hi_tmp[i_hi+1, j_hi+1] = lo[i_lo, j_lo]
            end
        end

        hi_tmp[1, :] = hi_tmp[end-1, :]
        hi_tmp[end, :] = hi_tmp[2, :]
        hi_tmp[:, 1] = hi_tmp[:, end-1]
        hi_tmp[:, end] = hi_tmp[:, 2]

        # hi[:, :] = hi_tmp[:, :]

        @tturbo for j in 2:size(hi, 2)-1
            for i in 2:size(hi, 1)-1
                hi[i, j] = ( 1//16*hi_tmp[i-1, j-1] + 1//8*hi_tmp[i, j-1] + 1//16*hi_tmp[i+1, j-1]
                           + 1// 8*hi_tmp[i-1, j  ] + 1//4*hi_tmp[i, j  ] + 1// 8*hi_tmp[i+1, j  ]
                           + 1//16*hi_tmp[i-1, j+1] + 1//8*hi_tmp[i, j+1] + 1//16*hi_tmp[i+1, j+1] )
            end
        end
    elseif n_hi ÷ n_lo == 3
        @tturbo for j in 0:size(lo, 2)-3
            for i in 0:size(lo, 1)-3
                i_lo = i + 2
                j_lo = j + 2
                i_hi = i*3 + 2
                j_hi = j*3 + 2
                hi_tmp[i_hi  , j_hi  ] = lo[i_lo, j_lo]
                hi_tmp[i_hi+1, j_hi  ] = lo[i_lo, j_lo]
                hi_tmp[i_hi+2, j_hi  ] = lo[i_lo, j_lo]
                hi_tmp[i_hi  , j_hi+1] = lo[i_lo, j_lo]
                hi_tmp[i_hi+1, j_hi+1] = lo[i_lo, j_lo]
                hi_tmp[i_hi+2, j_hi+1] = lo[i_lo, j_lo]
                hi_tmp[i_hi  , j_hi+2] = lo[i_lo, j_lo]
                hi_tmp[i_hi+1, j_hi+2] = lo[i_lo, j_lo]
                hi_tmp[i_hi+2, j_hi+2] = lo[i_lo, j_lo]
            end
        end

        hi_tmp[1, :] = hi_tmp[end-1, :]
        hi_tmp[end, :] = hi_tmp[2, :]
        hi_tmp[:, 1] = hi_tmp[:, end-1]
        hi_tmp[:, end] = hi_tmp[:, 2]

        @tturbo for j in 2:size(hi, 2)-1
            for i in 2:size(hi, 1)-1
                hi[i, j] = ( 1//9*hi_tmp[i-1, j-1] + 1//9*hi_tmp[i, j-1] + 1//9*hi_tmp[i+1, j-1]
                           + 1//9*hi_tmp[i-1, j  ] + 1//9*hi_tmp[i, j  ] + 1//9*hi_tmp[i+1, j  ]
                           + 1//9*hi_tmp[i-1, j+1] + 1//9*hi_tmp[i, j+1] + 1//9*hi_tmp[i+1, j+1] )
            end
        end
    else
        throw(DomainError(n_hi/n_lo, "Refinement should be 2 or 3"))
    end
end


## Set up the grids.
n_hi = 36
n_lo = n_hi ÷ 3

a_lo = rand(n_lo, n_lo)
a_hi = zeros(n_hi, n_hi)

a_hi_gc = zeros(n_hi+2, n_hi+2)
a_lo_gc = zeros(n_lo+2, n_lo+2)

a_hi_gc[2:n_hi+1, 2:n_hi+1] = a_hi[:, :]
a_lo_gc[2:n_lo+1, 2:n_lo+1] = a_lo[:, :]

a_hi_gc[1, :] = a_hi_gc[end-1, :]
a_hi_gc[end, :] = a_hi_gc[2, :]
a_hi_gc[:, 1] = a_hi_gc[:, end-1]
a_hi_gc[:, end] = a_hi_gc[:, 2]

a_lo_gc[1, :] = a_lo_gc[end-1, :]
a_lo_gc[end, :] = a_lo_gc[2, :]
a_lo_gc[:, 1] = a_lo_gc[:, end-1]
a_lo_gc[:, end] = a_lo_gc[:, 2]

dx_hi = 1/n_hi; dx_lo = 1/n_lo
x_hi = dx_hi/2:dx_hi:1 |> collect
x_lo = dx_lo/2:dx_lo:1 |> collect


## Compute.
a_hi_tmp = zeros(size(a_hi_gc))
@btime refine_field_int!(a_hi_gc, a_hi_tmp, a_lo_gc, n_hi, n_lo)
a_hi_int = a_hi_gc[2:end-1, 2:end-1]

println(mean(a_lo) ≈ mean(a_hi_int))


## Plot the output.
figure()
plot(x_hi, a_hi_int[1, :], "C0-o")
plot(x_hi, a_hi_int[2, :], "C1-+")
if n_hi ÷ n_lo == 3
    plot(x_hi, a_hi_int[3, :], "C2-^")
end
plot(x_lo, a_lo[1, :], "k:")
display(gcf())

xh_hi = 0:dx_hi:1 |> collect
xh_lo = 0:dx_lo:1 |> collect

figure(figsize=(10, 5))
subplot(121)
pcolormesh(xh_lo, xh_lo, a_lo, vmin=0, vmax=1)
subplot(122)
pcolormesh(xh_hi, xh_hi, a_hi_int, vmin=0, vmax=1)
tight_layout()
display(gcf())
show()
