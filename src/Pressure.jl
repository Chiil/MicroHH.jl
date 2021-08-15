using FFTW
using Printf
using HDF5

struct Pressure
    fft_forward 
    fft_backward
    bmati
    bmatj
    a
    c
end

function Pressure(g::Grid)
    nthreads = (Threads.nthreads() == 1) ? 1 : 2*Threads.nthreads()
    FFTW.set_num_threads(nthreads)

    tmp = rand(g.itot, g.jtot, g.ktot)
    fft_plan_f = FFTW.plan_r2r(tmp, FFTW.R2HC, [1, 2], flags=FFTW.MEASURE)
    fft_plan_b = FFTW.plan_r2r(tmp, FFTW.HC2R, [1, 2], flags=FFTW.MEASURE)

    bmati = zeros(g.itot)
    bmatj = zeros(g.jtot)
    a = zeros(g.ktot)
    c = zeros(g.ktot)

    work2d = zeros(g.itot, g.jtot)

    dxidxi = g.dxi^2
    dyidyi = g.dyi^2

    for j in 0:g.jtot÷2
        bmatj[j+1] = 2. * (cos(2pi*j/g.jtot) - 1.) * dyidyi;
    end

    for j in g.jtot÷2+1:g.jtot-1
        bmatj[j+1] = bmatj[g.jtot-j+1];
    end

    for i in 0:g.itot÷2
        bmati[i+1] = 2. * (cos(2pi*i/g.itot) - 1.) * dxidxi;
    end

    for i in g.itot÷2+1:g.itot-1
        bmati[i+1] = bmati[g.itot-i+1];
    end

    for k in 1:g.ktot
        a[k] = g.dz[k+g.kgc] * g.dzhi[k+g.kgc  ];
        c[k] = g.dz[k+g.kgc] * g.dzhi[k+g.kgc+1];
    end

    Pressure(fft_plan_f, fft_plan_b, bmati, bmatj, a, c)
end

function divergence_kernel(
    u, v, w,
    dxi, dyi, dzi,
    is, ie, js, je, ks, ke)

    divmax = 0.
    @inbounds for k in ks:ke
        for j in js:je
            for i in is:ie
                div = @fd (u, v, w) gradx(u) + grady(v) + gradz(w)
                divmax = max(divmax, div)
            end
        end
    end

    return divmax
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
    b[:, :, 1] .+= a[1]
    b[:, :, ktot] .+= c[ktot]

    # Set the wave number 0 to Dirichlet = 0.
    b[1, 1, ktot] -= 2c[ktot]
end

function solve_tdma_kernel!(
    p, work3d, work2d,
    a, b, c,
    ktot)

    @inbounds @. work2d[:, :] = b[:, :, 1]
    @inbounds @. p[:, :, 1] /= work2d[:, :]

    for k in 2:ktot
        @inbounds @. work3d[:, :, k] = c[k-1] / work2d[:, :]
        @inbounds @. work2d[:, :] = b[:, :, k] - a[k]*work3d[:, :, k]
        @inbounds @. p[:, :, k] -= a[k] * p[:, :, k-1]
        @inbounds @. p[:, :, k] /= work2d[:, :]
    end

    for k in ktot-1:1:-1
        @inbounds @. p[:, :, k] -= work3d[:, :, k+1] * p[:, :, k+1]
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
    f.w_tend[:, :, g.ks ] .= 0.
    f.w_tend[:, :, g.keh] .= 0.

    p_nogc = zeros(g.itot, g.jtot, g.ktot)

    dti_sub = 1/get_sub_dt(t)

    input_kernel!(
        p_nogc,
        f.u, f.v, f.w,
        f.u_tend, f.v_tend, f.w_tend,
        g.dxi, g.dyi, g.dzi, dti_sub,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    p_fft = p.fft_forward * p_nogc

    b = zeros(g.itot, g.jtot, g.ktot)
    solve_pre_kernel!(
        p_fft, b,
        g.dz, p.bmati, p.bmatj, p.a, p.c,
        g.itot, g.jtot, g.ktot, g.kgc)

    work3d = zeros(g.itot, g.jtot, g.ktot)
    work2d = zeros(g.itot, g.jtot)

    solve_tdma_kernel!(
        p_fft, work3d, work2d,
        p.a, b, p.c,
        g.ktot)

    p_nogc = (p.fft_backward * p_fft) ./ (g.itot * g.jtot)

    solve_post_kernel!(
        f.p, p_nogc,
        g.is, g.ie, g.js, g.je, g.ks, g.ke,
        g.igc, g.jgc, g.kgc)

    # Set the bot and top bc's
    @inbounds @. f.p[:, :, g.ks-1] = f.p[:, :, g.ks]
    @inbounds @. f.p[:, :, g.ke+1] = f.p[:, :, g.ke]

    boundary_cyclic_kernel!(
        f.p, g.is, g.ie, g.js, g.je, g.igc, g.jgc)

    output_kernel!(
        f.u_tend, f.v_tend, f.w_tend,
        f.p,
        g.dxi, g.dyi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end

function calc_divergence(f::Fields, g::Grid)
    div = divergence_kernel(
        f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
