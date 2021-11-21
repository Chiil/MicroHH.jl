## Loading packages.
using MicroHH
using PyPlot
using MicroHH.StencilBuilder
using LoopVectorization


## Loading settings.
include("drycbl_settings_post.jl")


## Initialize the model.
n_domains = 2
m = Model("drycbl", n_domains, settings, float_type)
load_model!(m)
in_progress = prepare_model!(m)


## Compute the scalar dissipation
function calc_scalar_dissipation!(
    slngrad, s,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (slngrad, s) slngrad = log(gradx(gradx(s))^2 + grady(grady(s))^2 + gradz(gradz(s))^2)
    end
end

f2 = m.fields[2]; g2 = m.grid[2]
slngrad = zeros(size(f2.s))
calc_scalar_dissipation!(
    slngrad, f2.s,
    g2.dxi, g2.dyi, g2.dzi, g2.dzhi,
    g2.is, g2.ie, g2.js, g2.je, g2.ks, g2.ke)

x = @view g2.x[g2.is:g2.ie]
y = @view g2.y[g2.js:g2.je]
z = @view g2.z[g2.ks:g2.ke]
slngrad_xy = @view slngrad[g2.is:g2.ie, g2.js:g2.je, 50]
slngrad_xz = @view slngrad[g2.is:g2.ie, g2.js, g2.ks:g2.ke]

smin = -25; smax = -11

figure()
pcolormesh(x, y, slngrad_xy', vmin=smin, vmax=smax, shading="nearest", cmap=plt.cm.magma)
xlabel("x (m)")
ylabel("y (m)")
colorbar()
tight_layout()
display(gcf())

figure()
pcolormesh(x, z, slngrad_xz', vmin=smin, vmax=smax, shading="nearest", cmap=plt.cm.magma)
xlabel("x (m)")
ylabel("z (m)")
colorbar()
tight_layout()
display(gcf())

show()