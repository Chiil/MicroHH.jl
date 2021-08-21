using FFTW
using Printf
using HDF5

struct Pressure{TF <: Union{Float32, Float64}}
    fft_forward
    fft_backward
    bmati::Vector{TF}
    bmatj::Vector{TF}
    a::Vector{TF}
    c::Vector{TF}
    work3d::Array{TF, 3}
    work2d::Array{TF, 2}
    b::Array{TF, 3}
    p_nogc
end

function Pressure(g::Grid, TF)
    nthreads = (Threads.nthreads() == 1) ? 1 : 2*Threads.nthreads()
    FFTW.set_num_threads(nthreads)

    # We set the type of rand here explictly to trigger the right precision in FFTW
    tmp = rand(TF, g.itot, g.jtot, g.ktot)
    fft_plan_f = FFTW.plan_r2r(tmp, FFTW.R2HC, [1, 2], flags=FFTW.MEASURE)
    fft_plan_b = FFTW.plan_r2r(tmp, FFTW.HC2R, [1, 2], flags=FFTW.MEASURE)

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

    work3d = zeros(g.itot, g.jtot, g.ktot)
    work2d = zeros(g.itot, g.jtot)
    b = zeros(g.itot, g.jtot, g.ktot)
    p_nogc = zeros(g.itot, g.jtot, g.ktot)

    Pressure{TF}(fft_plan_f, fft_plan_b, bmati, bmatj, a, c, work3d, work2d, b, p_nogc)
end

function input_kernel!(
    p,
    u, v, w,
    ut, vt, wt,
    dxi, dyi, dzi, dti,
    is, ie, js, je, ks, ke)

    @fast3d begin
        p[i-is+1, j-js+1, k-ks+1] = @fd (u, v, w, ut, vt, wt) gradx(ut + u*dti) + grady(vt + v*dti) + gradz(wt + w*dti)
    end
end

function solve_pre_kernel!(
    p, b,
    dz, bmati, bmatj, a, c,
    itot, jtot, ktot, kgc)

    @tturbo unroll=8 for k in 1:ktot
        for j in 1:jtot
            for i in 1:itot
                b[i, j, k] = dz[k+kgc]^2 * (bmati[i]+bmatj[j]) - (a[k]+c[k]);
                p[i, j, k] *= dz[k+kgc]^2
            end
        end
    end

    # Set the BCs of all wave numbers to Neumann = 0.
    @tturbo for j in 1:jtot
        for i in 1:itot
            b[i, j, 1] += a[1]
            b[i, j, ktot] += c[ktot]
        end
    end

    # Set the wave number 0 to Dirichlet = 0.
    b[1, 1, ktot] -= 2c[ktot]
end

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
        @fd (ut, p) ut -= gradx(p)
        @fd (vt, p) vt -= grady(p)
        @fd (wt, p) wt -= gradz(p)
    end
end

function calc_pressure_tend!(f::Fields, g::Grid, t::Timeloop, p::Pressure)
    # Set the cyclic boundaries for the tendencies.
    boundary_cyclic_kernel!(
        f.u_tend, g.is, g.ie, g.js, g.je, g.igc, g.jgc)
    boundary_cyclic_kernel!(
        f.v_tend, g.is, g.ie, g.js, g.je, g.igc, g.jgc)

    # Set the boundaries of wtend to zero.
    # CvH: Fix this later.
    @tturbo for j in 1:g.jcells
        for i in 1:g.icells
            f.w_tend[i, j, g.ks ] = 0
            f.w_tend[i, j, g.keh] = 0
        end
    end


    dti_sub = 1/get_sub_dt(t)

    input_kernel!(
        p.p_nogc,
        f.u, f.v, f.w,
        f.u_tend, f.v_tend, f.w_tend,
        g.dxi, g.dyi, g.dzi, dti_sub,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    p_fft = p.fft_forward * p.p_nogc

    solve_pre_kernel!(
        p_fft, p.b,
        g.dz, p.bmati, p.bmatj, p.a, p.c,
        g.itot, g.jtot, g.ktot, g.kgc)

    solve_tdma_kernel!(
        p_fft, p.work3d, p.work2d,
        p.a, p.b, p.c,
        g.itot, g.jtot, g.ktot)

    p.p_nogc[:, :, :] = (p.fft_backward * p_fft) ./ (g.itot * g.jtot)

    solve_post_kernel!(
        f.p, p.p_nogc,
        g.is, g.ie, g.js, g.je, g.ks, g.ke,
        g.igc, g.jgc, g.kgc)

    # Set the bot and top bc's
    @tturbo for j in 1:g.jcells
        for i in 1:g.icells
            f.p[i, j, g.ks-1] = f.p[i, j, g.ks]
            f.p[i, j, g.ke+1] = f.p[i, j, g.ke]
        end
    end

    boundary_cyclic_kernel!(
        f.p, g.is, g.ie, g.js, g.je, g.igc, g.jgc)

    output_kernel!(
        f.u_tend, f.v_tend, f.w_tend,
        f.p,
        g.dxi, g.dyi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end