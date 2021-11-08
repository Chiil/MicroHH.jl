using Interpolations
using BenchmarkTools

function interpolate!(a_hi, a_lo, x_hi, x_lo)
    interp_a = interpolate((x_lo, x_lo, x_lo), a_lo, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    a_hi[:, :, :] .= interp_a(x_hi, x_hi, x_hi)
end

n_hi = 256
n_lo = 128

dx_hi = 1//n_hi
dx_lo = 1//n_lo

x_hi = 1//2*dx_hi:dx_hi:1
x_lo = -1//2*dx_lo:dx_lo:1+1//2*dx_lo

a0_lo = rand(n_lo+2, n_lo+2, n_lo+2)
a0_hi = zeros(n_hi, n_hi, n_hi)

a1_lo = rand(n_lo+2, n_lo+2, n_lo+2)
a1_hi = zeros(n_hi, n_hi, n_hi)

a2_lo = rand(n_lo+2, n_lo+2, n_lo+2)
a2_hi = zeros(n_hi, n_hi, n_hi)

a3_lo = rand(n_lo+2, n_lo+2, n_lo+2)
a3_hi = zeros(n_hi, n_hi, n_hi)

@btime begin
    @sync begin
        Threads.@spawn interpolate!(a0_hi, a0_lo, x_hi, x_lo)
        Threads.@spawn interpolate!(a1_hi, a1_lo, x_hi, x_lo)
        Threads.@spawn interpolate!(a2_hi, a2_lo, x_hi, x_lo)
        Threads.@spawn interpolate!(a3_hi, a3_lo, x_hi, x_lo)
    end
end
