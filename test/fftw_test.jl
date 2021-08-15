using FFTW
using BenchmarkTools
using HDF5
using Tullio

itot = 32; jtot = 48; ktot = 64;

x = collect(0:itot-1) ./ itot .* 2pi
y = collect(0:jtot-1) ./ jtot .* 2pi
z = collect(0:ktot-1) ./ ktot .* 2pi

a = zeros(itot, jtot, ktot)
@tullio a[i, j, k] = cos(4x[i]) + cos(8y[j])

tmp = rand(itot, jtot, ktot)

FFTW.set_num_threads(Threads.nthreads())

fft_plan_f = FFTW.plan_r2r(a, FFTW.R2HC, (1, 2), flags=FFTW.ESTIMATE)
fft_plan_b = FFTW.plan_r2r(a, FFTW.HC2R, (1, 2), flags=FFTW.ESTIMATE)

tmp[:, :, :] = fft_plan_f * a
a_new = fft_plan_b * tmp ./ (itot * jtot)

h5open("fft_test.h5", "w") do fid
    write(fid, "tmp", tmp[:, :, :])
    write(fid, "a", a[:, :, :])
    write(fid, "a_new", a_new[:, :, :])
end
