## Packages.
using LoopVectorization
using BenchmarkTools

## Functions.
function solve_tdma_kernel!(
    p, work3d, work2d,
    a, b, c,
    itot, jtot, ktot)

    @tturbo for j in 1:jtot
        for i in 1:itot
            work2d[i, j] = b[i, j, 1]
            p[i, j, 1] /= work2d[i, j]
        end
    end

    for k in 2:ktot
        @tturbo for j in 1:jtot
            for i in 1:itot
                work3d[i, j, k] = c[k-1] / work2d[i, j]
                work2d[i, j] = b[i, j, k] - a[k]*work3d[i, j, k]
                p[i, j, k] -= a[k] * p[i, j, k-1]
                p[i, j, k] /= work2d[i, j]
            end
        end
    end

    for k in ktot-1:-1:1
        @tturbo for j in 1:jtot
            for i in 1:itot
                p[i, j, k] -= work3d[i, j, k+1] * p[i, j, k+1]
            end
        end
    end
end

function solve_tdma_kernel_2!(
    p, work3d, work2d,
    a, b, c,
    itot, jtot, ktot)

    @tturbo for j in 1:jtot
        for i in 1:itot
            work2d[i, j] = b[i, j, 1]
            p[i, j, 1] /= work2d[i, j]
        end
    end

    for k in 2:ktot
        @tturbo for j in 1:jtot
            for i in 1:itot
                work3d[i, j, k] = c[k-1] / work2d[i, j]
                work2d[i, j] = b[i, j, k] - a[k]*work3d[i, j, k]
                p[i, j, k] -= a[k] * p[i, j, k-1]
                p[i, j, k] /= work2d[i, j]
            end
        end
    end

    for k in ktot-1:-1:1
        @tturbo for j in 1:jtot
            for i in 1:itot
                p[i, j, k] -= work3d[i, j, k+1] * p[i, j, k+1]
            end
        end
    end
end

## Settings.
itot = 512; jtot = 1; ktot = 512;
p_ref = rand(itot, jtot, ktot) .+ 1
work3d = zeros(itot, jtot, ktot)
work2d = zeros(itot, jtot)
a = rand(ktot) .+ 1
b = rand(itot, jtot, ktot) .+ 1
c = rand(ktot) .+ 1
p = copy(p_ref)
p2 = copy(p_ref)

## First solver.
# @btime solve_tdma_kernel!($p, $work3d, $work2d, $a, $b, $c, $itot, $jtot, $ktot)

## Second solver.
@btime solve_tdma_kernel_2!($p2, $work3d, $work2d, $a, $b, $c, $itot, $jtot, $ktot)
