## Load the packages.
using Tullio
using LoopVectorization
using Statistics
using BenchmarkTools
using PyPlot
using Interpolations


## Create the refine function.
function refine_field_nn!(hi, lo, n_hi, n_lo)
    # if n_hi ÷ n_lo == 2
    #     @tullio hi[2(i-1)+l] = lo[i] (l in 1:2)
    # elseif n_hi ÷ n_lo == 3
    #     @tullio hi[3(i-1)+l] = lo[i] (l in 1:3)
    # else
    #     throw(DomainError(n_hi/n_lo, "Refinement should be 2 or 3"))
    # end

    if n_hi ÷ n_lo == 2
        @tturbo for i in 0:2:size(hi, 1)-3
            i_lo = i÷2+2
            i_hi = i+2
            hi[i_hi  ] = lo[i_lo]
            hi[i_hi+1] = lo[i_lo]
        end
    elseif n_hi ÷ n_lo == 3
        @tturbo for i in 0:3:size(hi, 1)-3
            i_lo = i÷3+2
            i_hi = i+2
            hi[i_hi  ] = lo[i_lo]
            hi[i_hi+1] = lo[i_lo]
            hi[i_hi+2] = lo[i_lo]
        end
    else
        throw(DomainError(n_hi/n_lo, "Refinement should be 2 or 3"))
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
    else
        throw(DomainError(n_hi/n_lo, "Refinement should be 2 or 3"))
    end
end

function refine_field_int2!(a_hi_gc, a_lo_gc, n_hi, n_lo)
    dx_hi = 1/n_hi; dx_lo = 1/n_lo
    x_hi = -dx_hi/2:dx_hi:1+dx_hi |> collect
    x_lo = -dx_lo/2:dx_lo:1+dx_lo |> collect
    interp = LinearInterpolation(x_lo, a_lo_gc)
    a_hi_gc[2:end-1] = interp(x_hi[2:end-1])
end

function refine_field_int3!(hi, hi_tmp, lo, n_hi, n_lo)
    if n_hi ÷ n_lo == 2
        @tturbo for i in 0:2:size(hi, 1)-3
            i_lo = i÷2+2
            i_hi = i+2
            hi_tmp[i_hi  ] = lo[i_lo]
            hi_tmp[i_hi+1] = lo[i_lo]
        end

        hi_tmp[1] = hi_tmp[end-1]
        hi_tmp[end] = hi_tmp[2]

        @tturbo for i in 2:size(hi, 1)-1
            hi[i] = 1//4*hi_tmp[i-1] + 1//2*hi_tmp[i] + 1//4*hi_tmp[i+1]
        end
    elseif n_hi ÷ n_lo == 3
        @tturbo for i in 0:3:size(hi, 1)-3
            i_lo = i÷3+2
            i_hi = i+2
            hi[i_hi  ] = lo[i_lo]
            hi[i_hi+1] = lo[i_lo]
            hi[i_hi+2] = lo[i_lo]
        end

        hi[1] = hi[end-1]
        hi[end] = hi[2]

        @tturbo for i in 2:size(hi, 1)-1
            hi[i] = 1//4*hi[i-1] + 1//2*hi[i] + 1//4*hi[i+1]
        end
    else
        throw(DomainError(n_hi/n_lo, "Refinement should be 2 or 3"))
    end
end


## Set up the grids.
n_hi = 32
n_lo = n_hi ÷ 2

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

dx_hi = 1/n_hi; dx_lo = 1/n_lo
x_hi = dx_hi/2:dx_hi:1 |> collect
x_lo = dx_lo/2:dx_lo:1 |> collect


## Compute.
@btime refine_field_nn!(a_hi_gc, a_lo_gc, n_hi, n_lo)
a_hi_nn = a_hi_gc[2:end-1]

@btime refine_field_int!(a_hi_gc, a_lo_gc, n_hi, n_lo)
a_hi_int = a_hi_gc[2:end-1]

@btime refine_field_int2!(a_hi_gc, a_lo_gc, n_hi, n_lo)
a_hi_int2 = a_hi_gc[2:end-1]

a_hi_tmp = zeros(size(a_hi_gc))
@btime refine_field_int3!(a_hi_gc, a_hi_tmp, a_lo_gc, n_hi, n_lo)
a_hi_int3 = a_hi_gc[2:end-1]

println(mean(a_lo) ≈ mean(a_hi_nn))
println(mean(a_lo) ≈ mean(a_hi_int))
println(mean(a_lo) ≈ mean(a_hi_int2))
println(mean(a_lo) ≈ mean(a_hi_int3))

println(mean(a_lo), " ", mean(a_hi_int), " ", mean(a_hi_int3))
println(a_hi_int3 .- a_hi_int)


## Plot the output.
figure()
#plot(x_hi, a_hi_nn, "C0-o")
#plot(x_hi, a_hi_int, "C1-o")
#plot(x_hi, a_hi_int2, "C2-^")
plot(x_hi, a_hi_int3, "C3-+")
plot(x_lo, a_lo, "k:")
display(gcf())
show()
