## Loading packages.
using MicroHH
using PyPlot
using MicroHH.StencilBuilder
using LoopVectorization


## Loading settings.
include("drycbl_settings.jl")
settings[1]["timeloop"]["start_time"] = settings[1]["timeloop"]["end_time"]


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings, float_type)
load_model!(m)


## Compute the scalar dissipation
function calc_scalar_dissipation!(
    slngrad, s,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (slngrad, s) slngrad = log(gradx(gradx(s))^2 + grady(grady(s))^2 + gradz(gradz(s))^2)
    end
end

f = m.fields[1]; g = m.grid[1]
slngrad = zeros(size(f.s))
calc_scalar_dissipation!(
    slngrad, f.s,
    g.dxi, g.dyi, g.dzi, g.dzhi,
    g.is, g.ie, g.js, g.je, g.ks, g.ke)

x = @view g.x[g.is:g.ie]
y = @view g.y[g.js:g.je]
z = @view g.z[g.ks:g.ke]
slngrad_xy = @view slngrad[g.is:g.ie, g.js:g.je, 6]
slngrad_xz = @view slngrad[g.is:g.ie, g.js, g.ks:g.ke]

smin = -25; smax = -11

figure()
#pcolormesh(x, y, slngrad_xy', vmin=smin, vmax=smax, shading="nearest", cmap=plt.cm.magma)
pcolormesh(x, y, slngrad_xy', shading="nearest", cmap=plt.cm.magma)
xlabel("x (m)")
ylabel("y (m)")
colorbar()
tight_layout()
display(gcf())

figure()
#pcolormesh(x, z, slngrad_xz', vmin=smin, vmax=smax, shading="nearest", cmap=plt.cm.magma)
pcolormesh(x, z, slngrad_xz', shading="nearest", cmap=plt.cm.magma)
xlabel("x (m)")
ylabel("z (m)")
colorbar()
tight_layout()
display(gcf())

show()
