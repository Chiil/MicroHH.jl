using LoopVectorization
using BenchmarkTools

function boundary_cyclic_kernel!(
    a, is, ie, js, je, igc, jgc)
    # East-west BCs
    @inbounds @. a[1:igc, :, :] = a[ie-igc+1:ie, :, :]
    @inbounds @. a[ie+1:end, :, :] = a[is:is+igc-1, :, :]

    # North-south BCs
    @inbounds @. a[:, 1:jgc, :] = a[:, je-jgc+1:je, :]
    @inbounds @. a[:, je+1:end, :] = a[:, js:js+igc-1, :]
end

function boundary_cyclic_kernel_loop!(
    a, is, ie, js, je, igc, jgc)
    # East-west BCs
    @tturbo for k in 1:size(a, 3)
        for j in 1:size(a, 2)
            for i in 1:igc
                a[i, j, k] = a[ie-i+1, j, k]
                a[ie+i, j, k] = a[is+i-1, j, k]
            end
        end
    end

    # North-south BCs
    @tturbo for k in 1:size(a, 3)
        for j in 1:jgc
            for i in 1:size(a, 1)
                a[i, j, k] = a[i, je-j+1, k]
                a[i, je+j, k] = a[i, js+j-1, k]
            end
        end
    end
end

function boundary_cyclic_kernel_turbo!(
    a, is, ie, js, je, igc, jgc)
    # East-west BCs
    @turbo @. a[1:igc, :, :] = a[ie-igc+1:ie, :, :]
    @turbo @. a[ie+1:end, :, :] = a[is:is+igc-1, :, :]

    # North-south BCs
    @turbo @. a[:, 1:jgc, :] = a[:, je-jgc+1:je, :]
    @turbo @. a[:, je+1:end, :] = a[:, js:js+igc-1, :]
end

a = rand(258, 3, 258)
is = 2; ie = 257; js = 2; je = 2; igc = 1; jgc = 1;

@btime boundary_cyclic_kernel!($a, $is, $ie, $js, $je, $igc, $jgc)
@btime boundary_cyclic_kernel_loop!($a, $is, $ie, $js, $je, $igc, $jgc)
@btime boundary_cyclic_kernel_turbo!($a, $is, $ie, $js, $je, $igc, $jgc)
