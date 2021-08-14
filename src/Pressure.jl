using FFTW

struct Pressure
    fft_forward 
    fft_backward
end

function Pressure(g::Grid)
    nthreads = (Threads.nthreads() == 1) ? 1 : 2*Threads.nthreads()
    FFTW.set_num_threads(nthreads)

    a = rand(g.itot, g.jtot, g.ktot)
    fft_plan_f = FFTW.plan_r2r(a, FFTW.R2HC, (1, 2), flags=FFTW.PATIENT)
    fft_plan_b = FFTW.plan_r2r(a, FFTW.HC2R, (1, 2), flags=FFTW.PATIENT)

    Pressure(fft_plan_f, fft_plan_b)
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
