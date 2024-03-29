using FFTW


struct Pressure{TF <: Union{Float32, Float64}}
    fft_forward_i
    fft_backward_i
    fft_forward_j
    fft_backward_j
    bmati::Vector{TF}
    bmatj::Vector{TF}
    a::Vector{TF}
    c::Vector{TF}
    work3d::Array{TF, 3}
    work2d::Array{TF, 2}
    b::Array{TF, 3}
    p_nogc::Array{TF, 3}
    sendbuf::Array{TF, 1}
    recvbuf::Array{TF, 1}
end


function Pressure(g::Grid, pp::Parallel, TF)
    nthreads = Threads.nthreads()
    FFTW.set_num_threads(nthreads)

    # We set the type of rand here explictly to trigger the right precision in FFTW
    tmp = rand(TF, g.itot, g.jmax, g.kblock)
    fft_plan_fi = FFTW.plan_r2r!(tmp, FFTW.R2HC, 1, flags=FFTW.MEASURE)
    fft_plan_bi = FFTW.plan_r2r!(tmp, FFTW.HC2R, 1, flags=FFTW.MEASURE)

    tmp = rand(TF, g.iblock, g.jtot, g.kblock)
    fft_plan_fj = FFTW.plan_r2r!(tmp, FFTW.R2HC, 2, flags=FFTW.MEASURE)
    fft_plan_bj = FFTW.plan_r2r!(tmp, FFTW.HC2R, 2, flags=FFTW.MEASURE)

    bmati = zeros(g.itot)
    bmatj = zeros(g.jtot)
    a = zeros(g.ktot)
    c = zeros(g.ktot)

    dxidxi = g.dxi^2
    dyidyi = g.dyi^2

    for j in 0:g.jtot÷2
        bmatj[j+1] = 2 * (cos(2pi*j/g.jtot) - 1) * dyidyi;
    end

    for j in g.jtot÷2+1:g.jtot-1
        bmatj[j+1] = bmatj[g.jtot-j+1];
    end

    for i in 0:g.itot÷2
        bmati[i+1] = 2 * (cos(2pi*i/g.itot) - 1) * dxidxi;
    end

    for i in g.itot÷2+1:g.itot-1
        bmati[i+1] = bmati[g.itot-i+1];
    end

    for k in 1:g.ktot
        a[k] = g.dz[k+g.kgc] * g.dzhi[k+g.kgc  ];
        c[k] = g.dz[k+g.kgc] * g.dzhi[k+g.kgc+1];
    end

    work3d = zeros(g.iblock, g.jblock, g.ktot)
    work2d = zeros(g.iblock, g.jblock)
    b = zeros(g.iblock, g.jblock, g.ktot)
    p_nogc = zeros(g.imax, g.jmax, g.ktot)

    if pp isa ParallelSerial
        sendbuf = zeros(0); recvbuf = zeros(0);
    elseif pp isa ParallelDistributed
        sendbuf = zeros(g.imax*g.jmax*g.ktot); recvbuf = zeros(g.imax*g.jmax*g.ktot);
    end

    Pressure{TF}(
        fft_plan_fi, fft_plan_bi, fft_plan_fj, fft_plan_bj,
        bmati, bmatj, a, c, work3d, work2d, b, p_nogc,
        sendbuf, recvbuf)
end


function input_kernel!(
    p,
    u, v, w,
    ut, vt, wt,
    dxi, dyi, dzi, dti,
    is, ie, js, je, ks, ke)

    @fast3d begin
        p[i-is+1, j-js+1, k-ks+1] = @fd (u[ih, j, k], v[i, jh, k], w[i, j, kh], ut[ih, j, k], vt[i, jh, k], wt[i, j, kh]) begin
            grad2x(ut + u*dti) + grad2y(vt + v*dti) + grad2z(wt + w*dti)
        end
    end
end


function solve_pre_kernel!(
    p, b,
    dz, bmati, bmatj, a, c,
    iblock, jblock, ktot, kgc,
    id_x, id_y)

    @tturbo unroll=8 for k in 1:ktot
        for j in 1:jblock
            for i in 1:iblock
                i_abs = id_y*iblock + i; j_abs = id_x*jblock + j
                b[i, j, k] = dz[k+kgc]^2 * (bmati[i_abs]+bmatj[j_abs]) - (a[k]+c[k]);
                p[i, j, k] *= dz[k+kgc]^2
            end
        end
    end

    # Set the BCs of all wave numbers to Neumann = 0.
    @tturbo for j in 1:jblock
        for i in 1:iblock
            b[i, j, 1] += a[1]
            b[i, j, ktot] += c[ktot]
        end
    end

    # Set the wave number 0 to Dirichlet = 0.
    if id_x == 0 && id_y == 0
        b[1, 1, ktot] -= 2*c[ktot]
    end
end


function solve_tdma_kernel!(
    p, work3d, work2d,
    a, b, c,
    iblock, jblock, ktot)

    @tturbo for j in 1:jblock
        for i in 1:iblock
            work2d[i, j] = b[i, j, 1]
            p[i, j, 1] /= work2d[i, j]
        end
    end

    for k in 2:ktot
        @tturbo for j in 1:jblock
            for i in 1:iblock
                work3d[i, j, k] = c[k-1] / work2d[i, j]
                work2d[i, j] = b[i, j, k] - a[k]*work3d[i, j, k]
                p[i, j, k] -= a[k] * p[i, j, k-1]
                p[i, j, k] /= work2d[i, j]
            end
        end
    end

    for k in ktot-1:-1:1
        @tturbo for j in 1:jblock
            for i in 1:iblock
                p[i, j, k] -= work3d[i, j, k+1] * p[i, j, k+1]
            end
        end
    end
end


function solve_post_kernel!(
    p, p_nogc,
    is, ie, js, je, ks, ke,
    igc, jgc, kgc)

    @fast3d begin
        p[i, j, k] = p_nogc[i-igc, j-jgc, k-kgc]
    end
end


function output_kernel!(
    ut, vt, wt,
    p,
    dxi, dyi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (ut[ih, j, k], p[i, j, k]) ut -= grad2x(p)
        @fd (vt[i, jh, k], p[i, j, k]) vt -= grad2y(p)
        @fd (wt[i, j, kh], p[i, j, k]) wt -= grad2z(p)
    end
end


function calc_pressure_tend!(
        f::Fields, g::Grid, t::Timeloop, p::Pressure, b::Boundary, pp::Parallel, to::TimerOutput)

    # Set the cyclic boundaries for the tendencies.
    @timeit to "boundary_cyclic_kernel" boundary_cyclic_kernel!(
        f.u_tend, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, pp)
    @timeit to "boundary_cyclic_kernel" boundary_cyclic_kernel!(
        f.v_tend, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, pp)

    # Set the boundaries of wtend to zero.
    # CvH: Fix this later.
    @timeit to "bot_top_bcs" @tturbo for j in 1:g.jcells
        for i in 1:g.icells
            f.w_tend[i, j, g.ks ] = 0
            f.w_tend[i, j, g.keh] = 0
        end
    end

    # Convert the dti to the type of the arrays to prevent casting in kernel.
    dti_sub = convert(typeof(f.u_tend[1]), 1/get_sub_dt(t))

    @timeit to "input_kernel" input_kernel!(
        p.p_nogc,
        f.u, f.v, f.w,
        f.u_tend, f.v_tend, f.w_tend,
        g.dxi, g.dyi, g.dzi, dti_sub,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    p_nogc_x = reshape(p.p_nogc, (g.itot, g.jmax, g.kblock))
    @timeit to "transpose_zx" transpose_zx!(p_nogc_x, p.p_nogc, p.sendbuf, p.recvbuf, g, pp)

    @timeit to "fft_forward_i" p_fft_tmp = p.fft_forward_i * p_nogc_x

    p_fft_tmp2 = reshape(p_fft_tmp, (g.iblock, g.jtot, g.kblock))
    @timeit to "transpose_xy" transpose_xy!(p_fft_tmp2, p_fft_tmp, p.sendbuf, p.recvbuf, g, pp)

    @timeit to "fft_forward_j" p_fft = p.fft_forward_j * p_fft_tmp2

    p_fft_z = reshape(p_fft, (g.iblock, g.jblock, g.ktot))
    @timeit to "transpose_yzt" transpose_yzt!(p_fft_z, p_fft, p.sendbuf, p.recvbuf, g, pp)

    @timeit to "solve_pre_kernel" solve_pre_kernel!(
        p_fft_z, p.b,
        g.dz, p.bmati, p.bmatj, p.a, p.c,
        g.iblock, g.jblock, g.ktot, g.kgc,
        pp.id_x, pp.id_y)

    @timeit to "solve_tdma_kernel" solve_tdma_kernel!(
        p_fft_z, p.work3d, p.work2d,
        p.a, p.b, p.c,
        g.iblock, g.jblock, g.ktot)

    @timeit to "transpose_zty" transpose_zty!(p_fft, p_fft_z, p.sendbuf, p.recvbuf, g, pp)

    @timeit to "fft_backward_j" p_fft_tmp2 = p.fft_backward_j * p_fft

    p_fft_tmp = reshape(p_fft_tmp2, (g.itot, g.jmax, g.kblock))
    @timeit to "transpose_yx" transpose_yx!(p_fft_tmp, p_fft_tmp2, p.sendbuf, p.recvbuf, g, pp)

    @timeit to "fft_backward_i" p_nogc_x = p.fft_backward_i * p_fft_tmp
    @tturbo p_nogc_x ./= (g.itot*g.jtot)

    @timeit to "transpose_xz" transpose_xz!(p.p_nogc, p_nogc_x, p.sendbuf, p.recvbuf, g, pp)

    @timeit to "solve_post_kernel" solve_post_kernel!(
        f.p, p.p_nogc,
        g.is, g.ie, g.js, g.je, g.ks, g.ke,
        g.igc, g.jgc, g.kgc)

    # Set the bot and top bc's
    @timeit to "bot_top_bcs" @tturbo for j in 1:g.jcells
        for i in 1:g.icells
            f.p[i, j, g.ks-1] = f.p[i, j, g.ks]
            f.p[i, j, g.ke+1] = f.p[i, j, g.ke]
        end
    end

    @timeit to "boundary_cyclic_kernel" boundary_cyclic_kernel!(
        f.p, g.is, g.ie, g.js, g.je, g.igc, g.jgc, b.buffers, pp)

    @timeit to "output_kernel" output_kernel!(
        f.u_tend, f.v_tend, f.w_tend,
        f.p,
        g.dxi, g.dyi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
