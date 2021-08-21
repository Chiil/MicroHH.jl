using LoopVectorization
using BenchmarkTools

function set_ghost_cells_bot_kernel!(
    a, a_bot, a_gradbot, dzh, ks)

    # fac = 0.5f0 # This works.
    fac = 1//2 # This does not work.

    @tturbo for j in 1:size(a, 2)
        for i in 1:size(a, 1)
            a[i, j, ks-1] = -a_gradbot[i, j]*dzh[ks] + a[i, j, ks]
            a_bot[i, j] = fac * (a[i, j, ks-1] + a[i, j, ks])
        end
    end
end

TF = Float32
# TF = Float64

a = rand(TF, 50, 50, 50)
a_bot = rand(TF, 50, 50)
a_gradbot = rand(TF, 50, 50)
dzh = rand(TF, 50)
ks = 2

@btime set_ghost_cells_bot_kernel!($a, $a_bot, $a_gradbot, $dzh, $ks)

