## Load the packages.
using PyPlot
using HDF5
using Statistics
using Interpolations


## Read the data.
x_01 = h5read("drycbl.01.00001800.h5", "x")
x_02 = h5read("drycbl.02.00001800.h5", "x")
y_01 = h5read("drycbl.01.00001800.h5", "y")
y_02 = h5read("drycbl.02.00001800.h5", "y")
z_01 = h5read("drycbl.01.00001800.h5", "z")
z_02 = h5read("drycbl.02.00001800.h5", "z")
zh_01 = h5read("drycbl.01.00001800.h5", "zh")
zh_02 = h5read("drycbl.02.00001800.h5", "zh")

s_bot_01 = h5read("drycbl.01.00001800.h5", "s_bot", (:, :))
s_bot_02 = h5read("drycbl.02.00001800.h5", "s_bot", (:, :))

s_bot_01 .-= mean(s_bot_01)
s_bot_02 .-= mean(s_bot_02)

println(zh_01[13], " ", zh_02[37])
w_xy_01 = h5read("drycbl.01.00001800.h5", "w", (:, :, 13)) # at 300 m
w_xy_02 = h5read("drycbl.02.00001800.h5", "w", (:, :, 37))

w_xz_01 = h5read("drycbl.01.00001800.h5", "w", (:, 1, :))
w_xz_02 = h5read("drycbl.02.00001800.h5", "w", (:, 1, :))


## Read for plotting.
xh_01 = h5read("drycbl.01.00001800.h5", "xh")
xh_02 = h5read("drycbl.02.00001800.h5", "xh")
yh_01 = h5read("drycbl.01.00001800.h5", "yh")
yh_02 = h5read("drycbl.02.00001800.h5", "yh")
push!(xh_01, 3200.)
push!(xh_02, 3200.)
push!(yh_01, 3200.)
push!(yh_02, 3200.)


## Surface temperature figure.
figure(figsize=(9,4.5))
subplot(121)
pcolormesh(xh_01, yh_01, s_bot_01', vmin=-0.6, vmax=0.6)
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(xh_02, yh_02, s_bot_02', vmin=-0.6, vmax=0.6)
xlabel("x (m)")
ylabel("y (m)")

tight_layout()
display(gcf())

println("variance in s_bot: ", var(s_bot_01), " ", var(s_bot_02))


## Surface temperature figure tweaked
s_01 = copy(s_bot_01)
s_02 = copy(s_bot_02)
@. s_01[s_01 < 0.1] = NaN
@. s_02[s_02 < 0.1] = NaN

figure(figsize=(5,4.5))
pcolormesh(xh_02, yh_02, s_02', vmin=0., vmax=0.6, cmap=PyPlot.cm.Reds, alpha=0.4)
pcolormesh(xh_01, yh_01, s_01', vmin=0., vmax=0.6, cmap=PyPlot.cm.Blues, alpha=0.4)
xlabel("x (m)")
ylabel("y (m)")

tight_layout()
display(gcf())



## Vertical velocity xy
figure(figsize=(9,4.5))
subplot(121)
pcolormesh(xh_01, yh_01, w_xy_01', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic)
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(xh_02, yh_02, w_xy_02', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic)
xlabel("x (m)")
ylabel("y (m)")

tight_layout()
display(gcf())

println("variance in w_xy: ", var(w_xy_01), " ", var(w_xy_02))


## Vertical velocity xy interpolated
interp = LinearInterpolation((x_02, y_02), w_xy_02, extrapolation_bc=Periodic())
w_xy_02_to_01 = interp(x_01, y_01)

"""
figure(figsize=(9,4.5))
subplot(121)
pcolormesh(xh_01, yh_01, w_xy_01', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic)
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(xh_01, yh_01, w_xy_02_to_01', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic)
xlabel("x (m)")
ylabel("y (m)")

tight_layout()
display(gcf())
"""


## Scatter xy
w_min = max(minimum(w_xy_01), minimum(w_xy_02_to_01))
w_max = min(maximum(w_xy_01), maximum(w_xy_02_to_01))

figure(figsize=(5.5, 4.5))
hexbin(vec(w_xy_01), vec(w_xy_02_to_01), gridsize=60, bins="log", cmap=plt.cm.cividis)
colorbar()
plot([-2, 5], [-2, 5], "w:")
xlim(w_min, w_max)
ylim(w_min, w_max)
xlabel(L"w_xy_01 (m s$^{-1}$)")
ylabel(L"w_xy_02 (m s$^{-1}$)")

tight_layout()
display(gcf())

println("correlation in w_xy: ", cor(vec(w_xy_01), vec(w_xy_02_to_01)))


## Vertical velocity xz
zp_01 = copy(z_01)
zp_02 = copy(z_02)
pushfirst!(zp_01, 0.)
pushfirst!(zp_02, 0.)
push!(zp_01, 3200.)
push!(zp_02, 3200.)

figure(figsize=(9,4.5))
subplot(121)
pcolormesh(xh_01, zp_01, w_xz_01', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic)
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(xh_02, zp_02, w_xz_02', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic)
xlabel("x (m)")
ylabel("y (m)")

tight_layout()
display(gcf())

show()
