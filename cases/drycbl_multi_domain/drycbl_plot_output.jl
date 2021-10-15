## Load the packages.
using PyPlot
using HDF5
using Statistics


## Read the data.
x_01 = h5read("drycbl.01.00005400.h5", "x")
x_02 = h5read("drycbl.02.00005400.h5", "x")
y_01 = h5read("drycbl.01.00005400.h5", "y")
y_02 = h5read("drycbl.02.00005400.h5", "y")
zh_01 = h5read("drycbl.01.00005400.h5", "zh")
zh_02 = h5read("drycbl.02.00005400.h5", "zh")

s_bot_01 = h5read("drycbl.01.00005400.h5", "s_bot", (:, :))
s_bot_02 = h5read("drycbl.02.00005400.h5", "s_bot", (:, :))

s_bot_01 .-= mean(s_bot_01)
s_bot_02 .-= mean(s_bot_02)

w_xy_01 = h5read("drycbl.01.00005400.h5", "w", (:, :, 13)) # at 300 m
w_xy_02 = h5read("drycbl.02.00005400.h5", "w", (:, :, 25))

w_xz_01 = h5read("drycbl.01.00005400.h5", "w", (:, 1, :))
w_xz_02 = h5read("drycbl.02.00005400.h5", "w", (:, 1, :))


## Surface temperature figure.
figure(figsize=(9,4.5))
subplot(121)
pcolormesh(x_01, y_01, s_bot_01, vmin=-0.5, vmax=0.5, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(x_02, y_02, s_bot_02, vmin=-0.5, vmax=0.5, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

tight_layout()

println(var(s_bot_01), " ", var(s_bot_02))


## Surface temperature figure tweaked
s_01 = copy(s_bot_01)
s_02 = copy(s_bot_02)
@. s_01[s_01 < 0.1] = NaN
@. s_02[s_02 < 0.1] = NaN

figure(figsize=(5,4.5))
pcolormesh(x_02, y_02, s_02, vmin=0., vmax=0.5, shading="nearest", cmap=PyPlot.cm.Reds, alpha=0.5)
pcolormesh(x_01, y_01, s_01, vmin=0., vmax=0.5, shading="nearest", cmap=PyPlot.cm.Blues, alpha=0.5)
xlabel("x (m)")
ylabel("y (m)")

tight_layout()



## Vertical velocity xy
figure(figsize=(9,4.5))
subplot(121)
pcolormesh(x_01, y_01, w_xy_01, vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(x_02, y_02, w_xy_02, vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

tight_layout()

println(var(w_xy_01), " ", var(w_xy_02))


## Vertical velocity xz
figure(figsize=(9,4.5))
subplot(121)
pcolormesh(x_01, zh_01, w_xz_01', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(x_02, zh_02, w_xz_02', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

tight_layout()


show()
display(gcf())
