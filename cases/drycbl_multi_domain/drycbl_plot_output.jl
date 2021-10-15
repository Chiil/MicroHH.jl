using PyPlot
using HDF5
using Statistics


## Read the data.
dx_01 = 3200 / 128
dx_02 = 3200 / 256
x_01 = dx_01/2:dx_01:3200 |> collect
x_02 = dx_02/2:dx_02:3200 |> collect
z_01 = 0:dx_01:3200 |> collect
z_02 = 0:dx_02:3200 |> collect

s_bot_01 = h5read("drycbl.01.00003600.h5", "s_bot", (:, :))
s_bot_02 = h5read("drycbl.02.00003600.h5", "s_bot", (:, :))

s_bot_01 .-= mean(s_bot_01)
s_bot_02 .-= mean(s_bot_02)

w_xy_01 = h5read("drycbl.01.00003600.h5", "w", (:, :, 13)) # at 300 m
w_xy_02 = h5read("drycbl.02.00003600.h5", "w", (:, :, 25))

w_xz_01 = h5read("drycbl.01.00003600.h5", "w", (:, 1, :))
w_xz_02 = h5read("drycbl.02.00003600.h5", "w", (:, 1, :))


## Surface temperature figure.
figure(figsize=(9,4.5))
subplot(121)
pcolormesh(x_01, x_01, s_bot_01, vmin=-0.5, vmax=0.5, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(x_02, x_02, s_bot_02, vmin=-0.5, vmax=0.5, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

tight_layout()

println(var(s_bot_01), " ", var(s_bot_02))


## Vertical velocity xy
figure(figsize=(9,4.5))
subplot(121)
pcolormesh(x_01, x_01, w_xy_01, vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(x_02, x_02, w_xy_02, vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

tight_layout()

println(var(w_xy_01), " ", var(w_xy_02))


## Vertical velocity xz
figure(figsize=(9,4.5))
subplot(121)
pcolormesh(x_01, z_01, w_xz_01', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

subplot(122)
pcolormesh(x_02, z_02, w_xz_02', vmin=-1.5, vmax=1.5, cmap=plt.cm.seismic, shading="nearest")
xlabel("x (m)")
ylabel("y (m)")

tight_layout()


show()
display(gcf())
