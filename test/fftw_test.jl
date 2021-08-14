using FFTW
using BenchmarkTools

itot = 256; jtot = 384; ktot = 192;

a = rand(itot, jtot, ktot)
tmp = rand(itot, jtot, ktot)

FFTW.set_num_threads(Threads.nthreads())

fft_plan_f = FFTW.plan_r2r(a, FFTW.R2HC, (1, 2), flags=FFTW.PATIENT)
fft_plan_b = FFTW.plan_r2r(a, FFTW.HC2R, (1, 2), flags=FFTW.PATIENT)

@btime begin
    tmp[:, :, :] = fft_plan_f * a
    a[:, :, :] = fft_plan_b * tmp
end
