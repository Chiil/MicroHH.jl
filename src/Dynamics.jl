## Packages
using LoopVectorization


## Stencil builder
function make_index(a, arrays, i, j, k)
    if a in arrays
        i += (a in [ Symbol("u"), Symbol("ut")]) ? 0.5 : 0
        j += (a in [ Symbol("v"), Symbol("vt")]) ? 0.5 : 0
        k += (a in [ Symbol("w"), Symbol("wt")]) ? 0.5 : 0

        if i < 0
            i_int = convert(Int, abs(i))
            ex_i = :( i-$i_int )
        elseif i > 0
            i_int = convert(Int, i)
            ex_i = :( i+$i_int )
        else
            ex_i = :( i )
        end

        if j < 0
            j_int = convert(Int, abs(j))
            ex_j = :( j-$j_int )
        elseif j > 0
            j_int = convert(Int, j)
            ex_j = :( j+$j_int )
        else
            ex_j = :( j )
        end

        if k < 0
            k_int = convert(Int, abs(k))
            ex_k = :( k-$k_int )
        elseif k > 0
            k_int = convert(Int, k)
            ex_k = :( k+$k_int )
        else
            ex_k = :( k )
        end

        return :( $a[ $ex_i, $ex_j, $ex_k] )
    else
        return :( $a )
    end
end


function process_expr(ex, arrays, i, j, k)
    ci1 = -1/16; ci2 =   9/16; ci3 =  9/16; ci4 = -1/16
    cg1 =  1/24; cg2 = -27/24; cg3 = 27/24; cg4 = -1/24

    n = 1

    if (isa(ex.args[1], Symbol) && ex.args[1] == Symbol("gradx"))
        ex.args[1] = Symbol("gradx_")
        ex = :( $ex * dxi )
    elseif (isa(ex.args[1], Symbol) && ex.args[1] == Symbol("grady"))
        ex.args[1] = Symbol("grady_")
        ex = :( $ex * dyi )
    elseif (isa(ex.args[1], Symbol) && ex.args[1] == Symbol("gradz"))
        ex.args[1] = Symbol("gradz_")
        ex = :( $ex * dzi )
    end

    args = ex.args
    while n <= length(args)
        if isa(args[n], Expr)
            args[n] = process_expr(args[n], arrays, i, j, k)
            n += 1
        elseif isa(args[n], Symbol)
            if args[n] == Symbol("gradx_")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i-1.5, j, k)
                    args[n+1] = process_expr(args[n+1], arrays, i-0.5, j, k)
                    args[n+2] = process_expr(args[n+2], arrays, i+0.5, j, k)
                    args[n+3] = process_expr(args[n+3], arrays, i+1.5, j, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i-1.5, j, k)
                    args[n+1] = make_index(args[n+1], arrays, i-0.5, j, k)
                    args[n+2] = make_index(args[n+2], arrays, i+0.5, j, k)
                    args[n+3] = make_index(args[n+3], arrays, i+1.5, j, k)
                end
                args[n  ] = :( $cg1 * $(args[n  ])  )
                args[n+1] = :( $cg2 * $(args[n+1])  )
                args[n+2] = :( $cg3 * $(args[n+2])  )
                args[n+3] = :( $cg4 * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("grady_")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i, j-1.5, k)
                    args[n+1] = process_expr(args[n+1], arrays, i, j-0.5, k)
                    args[n+2] = process_expr(args[n+2], arrays, i, j+0.5, k)
                    args[n+3] = process_expr(args[n+3], arrays, i, j+1.5, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i, j-1.5, k)
                    args[n+1] = make_index(args[n+1], arrays, i, j-0.5, k)
                    args[n+2] = make_index(args[n+2], arrays, i, j+0.5, k)
                    args[n+3] = make_index(args[n+3], arrays, i, j+1.5, k)
                end
                args[n  ] = :( $cg1 * $(args[n  ])  )
                args[n+1] = :( $cg2 * $(args[n+1])  )
                args[n+2] = :( $cg3 * $(args[n+2])  )
                args[n+3] = :( $cg4 * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("gradz_")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i, j, k-1.5)
                    args[n+1] = process_expr(args[n+1], arrays, i, j, k-0.5)
                    args[n+2] = process_expr(args[n+2], arrays, i, j, k+0.5)
                    args[n+3] = process_expr(args[n+3], arrays, i, j, k+1.5)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i, j, k-1.5)
                    args[n+1] = make_index(args[n+1], arrays, i, j, k-0.5)
                    args[n+2] = make_index(args[n+2], arrays, i, j, k+0.5)
                    args[n+3] = make_index(args[n+3], arrays, i, j, k+1.5)
                end
                args[n  ] = :( $cg1 * $(args[n  ])  )
                args[n+1] = :( $cg2 * $(args[n+1])  )
                args[n+2] = :( $cg3 * $(args[n+2])  )
                args[n+3] = :( $cg4 * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("interpx")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i-1.5, j, k)
                    args[n+1] = process_expr(args[n+1], arrays, i-0.5, j, k)
                    args[n+2] = process_expr(args[n+2], arrays, i+0.5, j, k)
                    args[n+3] = process_expr(args[n+3], arrays, i+1.5, j, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i-1.5, j, k)
                    args[n+1] = make_index(args[n+1], arrays, i-0.5, j, k)
                    args[n+2] = make_index(args[n+2], arrays, i+0.5, j, k)
                    args[n+3] = make_index(args[n+3], arrays, i+1.5, j, k)
                end
                args[n  ] = :( $ci1 * $(args[n  ])  )
                args[n+1] = :( $ci2 * $(args[n+1])  )
                args[n+2] = :( $ci3 * $(args[n+2])  )
                args[n+3] = :( $ci4 * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("interpy")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i, j-1.5, k)
                    args[n+1] = process_expr(args[n+1], arrays, i, j-0.5, k)
                    args[n+2] = process_expr(args[n+2], arrays, i, j+0.5, k)
                    args[n+3] = process_expr(args[n+3], arrays, i, j+1.5, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i, j-1.5, k)
                    args[n+1] = make_index(args[n+1], arrays, i, j-0.5, k)
                    args[n+2] = make_index(args[n+2], arrays, i, j+0.5, k)
                    args[n+3] = make_index(args[n+3], arrays, i, j+1.5, k)
                end
                args[n  ] = :( $ci1 * $(args[n  ])  )
                args[n+1] = :( $ci2 * $(args[n+1])  )
                args[n+2] = :( $ci3 * $(args[n+2])  )
                args[n+3] = :( $ci4 * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("interpz")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i, j, k-1.5)
                    args[n+1] = process_expr(args[n+1], arrays, i, j, k-0.5)
                    args[n+2] = process_expr(args[n+2], arrays, i, j, k+0.5)
                    args[n+3] = process_expr(args[n+3], arrays, i, j, k+1.5)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i, j, k-1.5)
                    args[n+1] = make_index(args[n+1], arrays, i, j, k-0.5)
                    args[n+2] = make_index(args[n+2], arrays, i, j, k+0.5)
                    args[n+3] = make_index(args[n+3], arrays, i, j, k+1.5)
                end
                args[n  ] = :( $ci1 * $(args[n  ])  )
                args[n+1] = :( $ci2 * $(args[n+1])  )
                args[n+2] = :( $ci3 * $(args[n+2])  )
                args[n+3] = :( $ci4 * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
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


macro fd(arrays, ex)
    i = (ex.args[1] in [ Symbol("u"), Symbol("ut")]) ? -0.5 : 0
    j = (ex.args[1] in [ Symbol("v"), Symbol("vt")]) ? -0.5 : 0
    k = (ex.args[1] in [ Symbol("w"), Symbol("wt")]) ? -0.5 : 0

    if isa(arrays, Symbol)
        arrays = :( [ $arrays ] )
    end

    ex = process_expr(ex, arrays.args, i, j, k)

    println("Generated stencil: ")
    println(ex)
    println("")

    return esc(ex)
end


macro fast3d(ex)
    ex_loop = quote
        @tturbo unroll=4 for k in ks:ke
            for j in js:je
                for i in is:ie
                    $ex
                end
            end
        end
    end
    return esc(ex_loop)
end


## Kernels.
function dynamics_u_kernel!(
    ut, u, v, w,
    visc,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (ut, u, v, w) ut += (
            - gradx(interpx(u) * interpx(u)) + visc * (gradx(gradx(u)))
            - grady(interpx(v) * interpy(u)) + visc * (grady(grady(u)))
            - gradz(interpx(w) * interpz(u)) + visc * (gradz(gradz(u))) )
    end
end


function dynamics_v_kernel!(
    vt, u, v, w,
    visc,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (vt, u, v, w) vt += (
            - gradx(interpy(u) * interpx(v)) + visc * (gradx(gradx(v)))
            - grady(interpy(v) * interpy(v)) + visc * (grady(grady(v)))
            - gradz(interpy(w) * interpz(v)) + visc * (gradz(gradz(v))) )
    end
end


function dynamics_w_kernel!(
    wt, u, v, w, s,
    visc,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (wt, u, v, w, s) wt += (
            - gradx(interpz(u) * interpx(w)) + visc * (gradx(gradx(w)))
            - grady(interpz(v) * interpy(w)) + visc * (grady(grady(w)))
            - gradz(interpz(w) * interpz(w)) + visc * (gradz(gradz(w)))
            + interpz(s) )
    end
end


function dynamics_s_kernel!(
    st, u, v, w, s,
    visc,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (st, u, v, w, s) st += (
            - gradx(u * interpx(s)) + visc * (gradx(gradx(s)))
            - grady(v * interpy(s)) + visc * (grady(grady(s)))
            - gradz(w * interpz(s)) + visc * (gradz(gradz(s))) )
    end
end


## Dynamics kernels launcher.
function calc_dynamics!(f::Fields, g::Grid)
    dynamics_u_kernel!(
        f.u_tend, f.u, f.v, f.w,
        f.visc,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    dynamics_v_kernel!(
        f.v_tend, f.u, f.v, f.w,
        f.visc,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    dynamics_w_kernel!(
        f.w_tend, f.u, f.v, f.w, f.s,
        f.visc,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    dynamics_s_kernel!(
        f.s_tend, f.u, f.v, f.w, f.s,
        f.visc,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
