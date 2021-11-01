## Load the packages.
using PyPlot
using HDF5
using Statistics
using Interpolations


## Read the data.
x_01 = h5read("drycbl.01.00003600.h5", "x")
x_02 = h5read("drycbl.02.00003600.h5", "x")
y_01 = h5read("drycbl.01.00003600.h5", "y")
y_02 = h5read("drycbl.02.00003600.h5", "y")
z_01 = h5read("drycbl.01.00003600.h5", "z")
z_02 = h5read("drycbl.02.00003600.h5", "z")
zh_01 = h5read("drycbl.01.00003600.h5", "zh")
zh_02 = h5read("drycbl.02.00003600.h5", "zh")

s_bot_01 = h5read("drycbl.01.00003600.h5", "s_bot", (:, :))
s_bot_02 = h5read("drycbl.02.00003600.h5", "s_bot", (:, :))

s_bot_01 .-= mean(s_bot_01)
s_bot_02 .-= mean(s_bot_02)

println(zh_01[17], " ", zh_02[49])
w_xy_01 = h5read("drycbl.01.00003600.h5", "w", (:, :, 25)) # at 300 m
w_xy_02 = h5read("drycbl.02.00003600.h5", "w", (:, :, 49))

w_xz_01 = h5read("drycbl.01.00003600.h5", "w", (:, 1, :))
w_xz_02 = h5read("drycbl.02.00003600.h5", "w", (:, 1, :))


## Read for plotting.
xh_01 = h5read("drycbl.01.00003600.h5", "xh")
xh_02 = h5read("drycbl.02.00003600.h5", "xh")
yh_01 = h5read("drycbl.01.00003600.h5", "yh")
yh_02 = h5read("drycbl.02.00003600.h5", "yh")
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

println(var(s_bot_01), " ", var(s_bot_02))


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

println(var(w_xy_01), " ", var(w_xy_02))


## Vertical velocity xy interpolated
interp = LinearInterpolation((x_02, y_02), w_xy_02, extrapolation_bc=Periodic())
w_xy_02_to_01 = interp(x_01, y_01)
println(size(w_xy_02_to_01))

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


## Vertical velocity xz
zp_01 = copy(z_01)
zp_02 = copy(z_02)
pushfirst!(zp_01, 0.)
pushfirst!(zp_02, 0.)
push!(zp_01, 3200.)
push!(zp_02, 3200.)

figure(figsize=(9,4.5))
subplot(121)
pcolormesh(xh_01, zh_01, w_xz_01', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic)
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(xh_02, zh_02, w_xz_02', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic)
xlabel("x (m)")
ylabel("y (m)")

tight_layout()
display(gcf())

show()
