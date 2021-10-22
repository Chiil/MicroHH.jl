## Load the packages.
using Tullio
using Statistics
using BenchmarkTools


## Set up the grids.
n_hi = 128
n_lo = 64

a_hi = rand(n_hi, n_hi, n_hi)
a_lo = zeros(n_lo, n_lo, n_lo)


## Coarsen the field
K = 1//8 * ones(2, 2, 2)
@tullio a_lo[i, j, k] = a_hi[2(i-1)+l, 2(j-1)+m, 2(k-1)+n] * K[l, m, n]


## Check the output.
println(mean(a_lo) ≈ mean(a_hi))
