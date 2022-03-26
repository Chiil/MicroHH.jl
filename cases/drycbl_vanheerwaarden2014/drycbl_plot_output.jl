## Load the packages.
using PyPlot
using HDF5
using Statistics
using Interpolations
using Printf
pygui(false)


## Read the data.
x = h5read("drycbl_cross.h5", "x")[:]
z = h5read("drycbl_cross.h5", "z")[:]
xh = [ h5read("drycbl_cross.h5", "xh")[:]; 2.]
zh = [ h5read("drycbl_cross.h5", "zh")[:]; 2.]
s = h5read("drycbl_cross.h5", "s_xz")

## Surface temperature figure.
for i in 50:20:size(s, 3)
    println("Processing frame $i")
    s_prime = s[:, :, i] .- mean(s[:, :, i], dims=1)

    fig = figure(figsize=(12.8, 4.0))
    pcolormesh(xh, zh, s_prime', vmin=-0.2, vmax=0.3, cmap=plt.cm.cividis)
    xlabel("x / δ")
    ylabel("z / δ")
    xlim(0., 2.0)
    ylim(0., 0.5)
    gca().set_aspect("equal")
    tight_layout()
    title("MicroHH.jl")

    filename = @sprintf("figs/fig_%05i.png", i)
    savefig(filename, dpi=100)
    close(fig)
end
