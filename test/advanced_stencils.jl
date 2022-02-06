## Packages.
using Tullio

include("../src/StencilBuilder.jl")
using .StencilBuilder


## Settings.
itot = 128; jtot = 96; ktot = 64
igc = 1; jgc = 1; kgc = 2
is = igc+1; js = jgc+1; ks = kgc+1
ie = igc+itot; je = jgc+jtot; ke = kgc+ktot

alpha = 1.
u = rand(itot+2igc, jtot+2jgc, ktot+2kgc)
v = rand(itot+2igc, jtot+2jgc, ktot+2kgc)
w = rand(itot+2igc, jtot+2jgc, ktot+2kgc)
s = rand(itot+2igc, jtot+2jgc, ktot+2kgc)
s_ref = rand(ktot+2kgc)

wt = rand(itot+2igc, jtot+2jgc, ktot+2kgc)

function parse_arrays(ex_arrays)
    arrays = Dict{Symbol, Tuple{Symbol, Symbol, Symbol}}()

    if !isa(ex_arrays, Expr) || !(ex_arrays.head in [ :vect, :tuple ])
        println("Array list $ex_arrays is invalid")
        return
    end

    for array in ex_arrays.args
        if !isa(array, Expr) || !(array.head == :ref)
            println("Array: $array is invalid")
            continue
        end

        array_name = array.args[1]
        array_dims = [:none, :none, :none]
        for arg in array.args[2:end]
            if arg == :i
                array_dims[1] = :ctr
            elseif arg == :ih
                array_dims[1] = :hlf
            elseif arg == :j
                array_dims[2] = :ctr
            elseif arg == :jh
                array_dims[2] = :hlf
            elseif arg == :k
                array_dims[3] = :ctr
            elseif arg == :kh
                array_dims[3] = :hlf
            else
                println("Array: $array is invalid")
            end
        end
        arrays[array_name] = tuple(array_dims...)
    end
    println("CvH final: $arrays")
end


## Settings.
# @fd_tullio (s, s_ref) (0, 0, -1//2) tmp = alpha * interpz(s - s_ref)
# @fd_tullio (wt[i, j, k], s[i, j, kh], s_ref[kh]) wt = alpha * interpz(s - s_ref)
parse_arrays(:( wt[i, j, k], s[i, j, kh], s_ref[kh], whker ))
