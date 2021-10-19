## Loading packages.
using MicroHH
using PyPlot
using MicroHH.StencilBuilder
using LoopVectorization


## Loading settings.
include("drycbl_settings_post.jl")


## Initialize the model.
n_domains = 2
m = Model("drycbl", n_domains, settings, Float32)
load_model!(m)
in_progress = prepare_model!(m)


## Compute the scalar dissipation
function calc_scalar_dissipation!(
    slngrad, s,
    dxi, dyi, dzi, dzhi,
    is, ie, js, je, ks, ke)

    @fast3d begin
        @fd (slngrad, s) s = gradx(gradx(s))^2 + grady(grady(s))^2 + gradz(gradz(s))^2
    end
end

f2 = m.fields[2]; g2 = m.grid[2]
slngrad = zeros(size(f2.s))
calc_scalar_dissipation!(
    slngrad, f2.s,
    g2.dxi, g2.dyi, g2.dzi, g2.dzhi,
    g2.is, g2.ie, g2.js, g2.je, g2.ks, g2.ke)
