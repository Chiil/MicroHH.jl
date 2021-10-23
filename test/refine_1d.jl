## Load the packages.
using Tullio
using Statistics
using BenchmarkTools
using PyPlot


## Create the refine function.
function refine_field!(hi, lo)
    @tullio hi[2(i-1)+l] = lo[i] (l in 1:2)
    return
end


## Set up the grids.
n_hi = 64
n_lo = 32

a_lo = rand(n_lo)
a_hi = zeros(n_hi)

a_hi_gc = zeros(n_hi+2)
a_lo_gc = zeros(n_lo+2)

a_hi_gc[2:n_hi+1] = a_hi[:]
a_lo_gc[2:n_lo+1] = a_lo[:]

a_hi_gc[1] = a_hi[end-1]
a_hi_gc[end] = a_hi[2]
a_lo_gc[1] = a_lo[end-1]
a_lo_gc[end] = a_lo[2]


## Compute.
a_hi_nn = copy(a_hi)
refine_field!(a_hi_nn, a_lo)
#refine_field!(a_hi_gc, a_lo_gc)

#a_hi[:] = a_hi_gc[2:end-1]
#a_lo[:] = a_lo_gc[2:end-1]


## Check the output.
dx_hi = 1/n_hi; dx_lo = 1/n_lo
x_hi = dx_hi/2:dx_hi:1 |> collect

dx_lo = 1/n_lo; dx_lo = 1/n_lo
x_lo = dx_lo/2:dx_lo:1 |> collect

figure()
plot(x_hi, a_hi_nn, "C0-o")
plot(x_lo, a_lo, "k:")
display(gcf())
show()

println(a_lo)
