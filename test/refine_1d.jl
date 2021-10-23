## Load the packages.
using Tullio
using LoopVectorization
using Statistics
using BenchmarkTools
using PyPlot


## Create the refine function.
function refine_field_nn!(hi, lo, n_hi, n_lo)
    if n_hi ÷ n_lo == 2
        @tullio hi[2(i-1)+l] = lo[i] (l in 1:2)
    elseif n_hi ÷ n_lo == 3
        @tullio hi[3(i-1)+l] = lo[i] (l in 1:3)
    end
end

function refine_field_int!(hi, lo, n_hi, n_lo)
    if n_hi ÷ n_lo == 2
        @tturbo for i in 0:2:size(hi, 1)-3
            i_lo = i÷2+2
            i_hi = i+2
            hi[i_hi  ] = 3//4*lo[i_lo] + 1//4*lo[i_lo-1]
            hi[i_hi+1] = 3//4*lo[i_lo] + 1//4*lo[i_lo+1]
        end
    elseif n_hi ÷ n_lo == 3
        @tturbo for i in 0:3:size(hi, 1)-3
            i_lo = i÷3+2
            i_hi = i+2
            hi[i_hi  ] = 2//3*lo[i_lo] + 1//3*lo[i_lo-1]
            hi[i_hi+1] = lo[i_lo]
            hi[i_hi+2] = 2//3*lo[i_lo] + 1//3*lo[i_lo+1]
        end
    end
end


## Set up the grids.
n_hi = 48
n_lo = 16

a_lo = rand(n_lo)
a_hi = zeros(n_hi)

a_hi_gc = zeros(n_hi+2)
a_lo_gc = zeros(n_lo+2)

a_hi_gc[2:n_hi+1] = a_hi[:]
a_lo_gc[2:n_lo+1] = a_lo[:]

a_hi_gc[1] = a_hi_gc[end-1]
a_hi_gc[end] = a_hi_gc[2]
a_lo_gc[1] = a_lo_gc[end-1]
a_lo_gc[end] = a_lo_gc[2]


## Compute.
a_hi_nn = copy(a_hi)
@btime refine_field_nn!(a_hi_nn, a_lo, n_hi, n_lo)

@btime refine_field_int!(a_hi_gc, a_lo_gc, n_hi, n_lo)
a_hi_int = a_hi_gc[2:end-1]


## Check the output.
dx_hi = 1/n_hi; dx_lo = 1/n_lo
x_hi = dx_hi/2:dx_hi:1 |> collect

dx_lo = 1/n_lo; dx_lo = 1/n_lo
x_lo = dx_lo/2:dx_lo:1 |> collect

figure()
plot(x_hi, a_hi_nn, "C0-o")
plot(x_hi, a_hi_int, "C1-o")
plot(x_lo, a_lo, "k:")
display(gcf())
show()

println(a_lo)
