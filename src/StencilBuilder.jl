module StencilBuilder


export @fast3d, @fd, @fd_tullio


function make_index(a, arrays, i, j, k)
    if haskey(arrays, a)
        if arrays[a][1] == :hlf i += 0.5 end
        if arrays[a][2] == :hlf j += 0.5 end
        if arrays[a][3] == :hlf k += 0.5 end

        ex_ijk = []

        if arrays[a][1] != :none
            if i < 0
                i_int = convert(Int, abs(i))
                ex_i = :( i-$i_int )
            elseif i > 0
                i_int = convert(Int, i)
                ex_i = :( i+$i_int )
            else
                ex_i = :( i+0 )
            end
            push!(ex_ijk, ex_i)
        end

        if arrays[a][2] != :none
            if j < 0
                j_int = convert(Int, abs(j))
                ex_j = :( j-$j_int )
            elseif j > 0
                j_int = convert(Int, j)
                ex_j = :( j+$j_int )
            else
                ex_j = :( j+0 )
            end
            push!(ex_ijk, ex_j)
        end

        if arrays[a][3] != :none
            if k < 0
                k_int = convert(Int, abs(k))
                ex_k = :( k-$k_int )
            elseif k > 0
                k_int = convert(Int, k)
                ex_k = :( k+$k_int )
            else
                ex_k = :( k+0 )
            end
            push!(ex_ijk, ex_k)
        end

        return :( $a[ $(ex_ijk...) ] )
    else
        return :( $a )
    end
end


function process_expr(ex, arrays, i, j, k)
    ci1 = 1//2; ci2 = 1//2;
    cg1 = -1; cg2 = 1;

    n = 1

    # Check whether expression is a gradient, and if so inject the appropriate dx/dy/dz.
    if (isa(ex.args[1], Symbol) && ex.args[1] == Symbol("gradx"))
        ex.args[1] = Symbol("gradx_")
        ex = :( $ex * dxi )
    elseif (isa(ex.args[1], Symbol) && ex.args[1] == Symbol("grady"))
        ex.args[1] = Symbol("grady_")
        ex = :( $ex * dyi )
    elseif (isa(ex.args[1], Symbol) && ex.args[1] == Symbol("gradz"))
        ex.args[1] = Symbol("gradz_")

        # Check whether the vertical location is at ctr or hlf location.
        # ctr location
        if isinteger(k)
            k_int = convert(Int, abs(k))
            if k > 0
                ex = :( $ex * dzi[k+$k_int] )
            elseif k < 0
                ex = :( $ex * dzi[k-$k_int] )
            else
                ex = :( $ex * dzi[k+0] )
            end
        # hlf location
        else
            k_int = convert(Int, abs(k + 1/2))
            if k > -1/2
                ex = :( $ex * dzhi[k+$k_int] )
            elseif k < -1/2
                ex = :( $ex * dzhi[k-$k_int] )
            else
                ex = :( $ex * dzhi[k+0] )
            end
        end
    end

    # Recurse through expression and replace gradients and interpolations by
    # the appropriate code.
    args = ex.args
    while n <= length(args)
        if isa(args[n], Expr)
            args[n] = process_expr(args[n], arrays, i, j, k)
            n += 1
        elseif isa(args[n], Symbol)
            if args[n] == Symbol("gradx_")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    args[n  ] = process_expr(args[n  ], arrays, i-0.5, j, k)
                    args[n+1] = process_expr(args[n+1], arrays, i+0.5, j, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    args[n  ] = make_index(args[n  ], arrays, i-0.5, j, k)
                    args[n+1] = make_index(args[n+1], arrays, i+0.5, j, k)
                end
                args[n  ] = :( $cg1 * $(args[n  ])  )
                args[n+1] = :( $cg2 * $(args[n+1])  )
                insert!(args, n, Symbol("+"))
                n += 3
            elseif args[n] == Symbol("grady_")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    args[n  ] = process_expr(args[n  ], arrays, i, j-0.5, k)
                    args[n+1] = process_expr(args[n+1], arrays, i, j+0.5, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    args[n  ] = make_index(args[n  ], arrays, i, j-0.5, k)
                    args[n+1] = make_index(args[n+1], arrays, i, j+0.5, k)
                end
                args[n  ] = :( $cg1 * $(args[n  ])  )
                args[n+1] = :( $cg2 * $(args[n+1])  )
                insert!(args, n, Symbol("+"))
                n += 3
            elseif args[n] == Symbol("gradz_")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    args[n  ] = process_expr(args[n  ], arrays, i, j, k-0.5)
                    args[n+1] = process_expr(args[n+1], arrays, i, j, k+0.5)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    args[n  ] = make_index(args[n  ], arrays, i, j, k-0.5)
                    args[n+1] = make_index(args[n+1], arrays, i, j, k+0.5)
                end
                args[n  ] = :( $cg1 * $(args[n  ])  )
                args[n+1] = :( $cg2 * $(args[n+1])  )
                insert!(args, n, Symbol("+"))
                n += 3
            elseif args[n] == Symbol("interpx")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    args[n  ] = process_expr(args[n  ], arrays, i-0.5, j, k)
                    args[n+1] = process_expr(args[n+1], arrays, i+0.5, j, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    args[n  ] = make_index(args[n  ], arrays, i-0.5, j, k)
                    args[n+1] = make_index(args[n+1], arrays, i+0.5, j, k)
                end
                args[n  ] = :( $ci1 * $(args[n  ])  )
                args[n+1] = :( $ci2 * $(args[n+1])  )
                insert!(args, n, Symbol("+"))
                n += 3
            elseif args[n] == Symbol("interpy")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    args[n  ] = process_expr(args[n  ], arrays, i, j-0.5, k)
                    args[n+1] = process_expr(args[n+1], arrays, i, j+0.5, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    args[n  ] = make_index(args[n  ], arrays, i, j-0.5, k)
                    args[n+1] = make_index(args[n+1], arrays, i, j+0.5, k)
                end
                args[n  ] = :( $ci1 * $(args[n  ])  )
                args[n+1] = :( $ci2 * $(args[n+1])  )
                insert!(args, n, Symbol("+"))
                n += 3
            elseif args[n] == Symbol("interpz")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    args[n  ] = process_expr(args[n  ], arrays, i, j, k-0.5)
                    args[n+1] = process_expr(args[n+1], arrays, i, j, k+0.5)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    args[n  ] = make_index(args[n  ], arrays, i, j, k-0.5)
                    args[n+1] = make_index(args[n+1], arrays, i, j, k+0.5)
                end
                args[n  ] = :( $ci1 * $(args[n  ])  )
                args[n+1] = :( $ci2 * $(args[n+1])  )
                insert!(args, n, Symbol("+"))
                n += 3
            else
                args[n] = make_index(args[n], arrays, i, j, k)
                n += 1
            end
        else
            n += 1
        end
    end
    return ex
end


function parse_arrays(ex_arrays)
    arrays = Dict{Symbol, Tuple{Symbol, Symbol, Symbol}}()

    if !isa(ex_arrays, Expr) || !(ex_arrays.head in [:vect, :tuple])
        @error "Array list $ex_arrays is invalid"
        return # CvH throw error
    end

    for array in ex_arrays.args
        if !isa(array, Expr) || !(array.head == :ref)
            @error "Array: $array is invalid"
            # CvH throw error
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
                @error "Array: $array is invalid"
                # CvH throw error
            end
        end
        arrays[array_name] = tuple(array_dims...)
    end
    return arrays
end


macro fd(ex_arrays, ex) build_fd(ex_arrays, :([0, 0, 0]), ex, false) end
macro fd(ex_arrays, ex_offset, ex) build_fd(ex_arrays, ex_offset, ex, false) end
macro fd_tullio(ex_arrays, ex) build_fd(ex_arrays, :([0, 0, 0]), ex, true) end
macro fd_tullio(ex_arrays, ex_offset, ex) build_fd(ex_arrays, ex_offset, ex, true) end


function build_fd(ex_arrays, ex_offset, ex, use_tullio)
    # Strip the begin end block from the expression in case present.
    if ex.head == :block
        ex = ex.args[2]
    end

    offset = eval(ex_offset)
    arrays = parse_arrays(ex_arrays)

    i = 0.; j = 0.; k = 0.
    if haskey(arrays, ex.args[1])
        if arrays[ex.args[1]][1] == :hlf i -= 0.5 end
        if arrays[ex.args[1]][2] == :hlf j -= 0.5 end
        if arrays[ex.args[1]][3] == :hlf k -= 0.5 end
    end

    i += offset[1]; j += offset[2]; k += offset[3]
        
    ex = process_expr(ex, arrays, i, j, k)

    if use_tullio
        ex = :( @tullio $ex i in is:ie, j in js:je, k in ks:ke )
    end

    @debug "Generated stencil:"
    @debug "$ex"
    @debug ""

    return esc(ex)
end


macro fast3d(ex)
    ex_loop = quote
        @tturbo unroll=8 for k in ks:ke
            for j in js:je
                for i in is:ie
                    $ex
                end
            end
        end
    end
    return esc(ex_loop)
end

end
