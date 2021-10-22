## Load the packages.
using Tullio
using Statistics
using BenchmarkTools


## Create the coarsening function.
function coarsen_field!(lo, hi)
    filter = 1//8 * ones(size(hi) .÷ size(lo))
    @tullio lo[i, j, k] = hi[2(i-1)+l, 2(j-1)+m, 2(k-1)+n] * filter[l, m, n]
    return
end


## Set up the grids.
n_hi = 128
n_lo = 64

a_hi = rand(n_hi, n_hi, n_hi)
a_lo = zeros(n_lo, n_lo, n_lo)


## Compute.
@btime coarsen_field!(a_lo, a_hi)


## Check the output.
println(mean(a_lo) ≈ mean(a_hi))
