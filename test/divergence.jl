## Packages
using BenchmarkTools
using LoopVectorization


## Reference kernel
function divergence_ref(
    u, v, w,
    is, ie, js, je, ks, ke)

    divmax = 0
    @inbounds for k in ks:ke
        for j in js:je
            for i in is:ie
                div = u[i, j, k] + v[i, j, k] + w[i, j, k]
                divmax = max(abs(div), divmax)
            end
        end
    end

    return divmax
end


## Reference kernel
function divergence(
    u, v, w,
    is, ie, js, je, ks, ke)

    divmax = 0
    @tturbo for k in ks:ke
        for j in js:je
            for i in is:ie
                div = u[i, j, k] + v[i, j, k] + w[i, j, k]
                divmax = max(abs(div), divmax)
            end
        end
    end

    return divmax
end


## User input
itot = 384; jtot = 384; ktot = 384;
igc = 1; jgc = 1; kgc = 1;

icells = itot + 2igc; jcells = jtot + 2jgc; kcells = ktot + 2kgc
is = igc+1; ie = igc+itot; js = jgc+1; je = jgc+jtot; ks = kgc+1; ke = kgc+ktot

u = rand(icells, jcells, kcells)
v = rand(icells, jcells, kcells)
w = rand(icells, jcells, kcells)


## Run reference kernel
div_ref = divergence_ref(
    u, v, w,
    is, ie, js, je, ks, ke)

@btime dummy = divergence_ref(
    $u, $v, $w,
    $is, $ie, $js, $je, $ks, $ke)


## Run fast kernel
div = divergence(
    u, v, w,
    is, ie, js, je, ks, ke)

@btime divergence(
    $u, $v, $w,
    $is, $ie, $js, $je, $ks, $ke)


## Compare results
println(div_ref, " ", div, " ", isapprox(div, div_ref))

