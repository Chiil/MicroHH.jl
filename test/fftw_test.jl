using FFTW
using BenchmarkTools
using HDF5
using Tullio

itot = 256; jtot = 192; ktot = 128;

float_type = Float32

x = collect(float_type, 0:itot-1) ./ itot .* 2pi
y = collect(float_type, 0:jtot-1) ./ jtot .* 2pi
z = collect(float_type, 0:ktot-1) ./ ktot .* 2pi

a = zeros(float_type, itot, jtot, ktot)
a_new = zeros(float_type, itot, jtot, ktot)
tmp = zeros(float_type, itot, jtot, ktot)

@tullio a[i, j, k] = cos(4x[i]) + cos(8y[j])


FFTW.set_num_threads(Threads.nthreads())

fft_plan_f = FFTW.plan_r2r(a, FFTW.R2HC, (1, 2), flags=FFTW.ESTIMATE)
fft_plan_b = FFTW.plan_r2r(a, FFTW.HC2R, (1, 2), flags=FFTW.ESTIMATE)

function forward_backward!(tmp, a_new, a, fft_plan_f, fft_plan_b, itot, jtot)
    tmp[:, :, :] = fft_plan_f * a
    a_new[:, :, :] = fft_plan_b * tmp ./ (itot * jtot)
end

@btime forward_backward!(
    $tmp, $a_new, $a, $fft_plan_f, $fft_plan_b, $itot, $jtot)
forward_backward!(tmp, a_new, a, fft_plan_f, fft_plan_b, itot, jtot)

h5open("fft_test.h5", "w") do fid
    write(fid, "tmp", tmp[:, :, :])
    write(fid, "a", a[:, :, :])
    write(fid, "a_new", a_new[:, :, :])
end
