## Loading packages.
using MicroHH
using PyPlot
using MicroHH.StencilBuilder
using LoopVectorization


## Loading settings.
include("drycbl_settings.jl")
for d in settings
    d["timeloop"]["start_time"] = 1800.
end


## Initialize the model.
n_domains = 2
m = Model("drycbl", n_domains, settings, float_type)
load_model!(m)
in_progress = prepare_model!(m)


## Compute the scalar dissipation.
function calc_scalar_dissipation!(
    slngrad, s,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (slngrad, s) slngrad = log(interpx(gradx(s))^2 + interpy(grady(s))^2 + interpz(gradz(s))^2)
    end
end

f2 = m.fields[2]; g2 = m.grid[2]
slngrad = zeros(size(f2.s))
calc_scalar_dissipation!(
    slngrad, f2.s,
    g2.dxi, g2.dyi, g2.dzi, g2.dzhi,
    g2.is, g2.ie, g2.js, g2.je, g2.ks, g2.ke)


## Plot the data.
x = @view g2.x[g2.is:g2.ie]
y = @view g2.y[g2.js:g2.je]
z = @view g2.z[g2.ks:g2.ke]
xh = @view g2.xh[g2.is:g2.ie+1]
yh = @view g2.yh[g2.js:g2.je+1]
slngrad_xy = @view slngrad[g2.is:g2.ie, g2.js:g2.je, 17]
slngrad_xz = @view slngrad[g2.is:g2.ie, g2.js, g2.ks:g2.ke]

smin = -15; smax = -9

figure()
pcolormesh(xh, yh, slngrad_xy', vmin=smin, vmax=smax, cmap=plt.cm.cividis)
xlabel("x (m)")
ylabel("y (m)")
colorbar()
tight_layout()
display(gcf())

figure()
pcolormesh(x, z, slngrad_xz', vmin=smin, vmax=smax, shading="nearest", cmap=plt.cm.cividis)
xlabel("x (m)")
ylabel("z (m)")
colorbar()
tight_layout()
display(gcf())

show()
