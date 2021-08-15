using FFTW

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
    fft_plan_f = FFTW.plan_r2r(tmp, FFTW.R2HC, (1, 2), flags=FFTW.PATIENT)
    fft_plan_b = FFTW.plan_r2r(tmp, FFTW.HC2R, (1, 2), flags=FFTW.PATIENT)

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

    println(bmati)
    println(bmatj)

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

function output_kernel!(
    ut, vt, wt,
    dxi, dyi, dzhi,
    is, ie, js, je, ks, ke)
end

function calc_pressure_tend!(f::Fields, g::Grid, t::Timeloop, p::Pressure)
    boundary_cyclic_kernel!(
        f.u_tend, g.is, g.ie, g.js, g.je, g.igc, g.jgc)
    boundary_cyclic_kernel!(
        f.v_tend, g.is, g.ie, g.js, g.je, g.igc, g.jgc)

    p_nogc = zeros(g.itot, g.jtot, g.ktot)

    input_kernel!(
        p_nogc,
        f.u, f.v, f.w,
        f.u_tend, f.v_tend, f.w_tend,
        g.dxi, g.dyi, g.dzi, 1/t.dt,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    tmp = p.fft_forward * p_nogc
    p_nogc = (p.fft_backward * tmp) ./ (g.itot * g.jtot)

    output_kernel!(
        f.u_tend, f.v_tend, f.w_tend,
        g.dxi, g.dyi, g.dzhi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end

function calc_divergence(f::Fields, g::Grid)
    div = divergence_kernel(
        f.u, f.v, f.w,
        g.dxi, g.dyi, g.dzi,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end
