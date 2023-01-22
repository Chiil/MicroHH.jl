## Loading packages.
using MicroHH
using Tullio
using Statistics


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings)


## Create the initials fields.
f = m.fields[1]; g = m.grid[1]
s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]
u = @view f.u[g.is:g.ie+1, g.js:g.je, g.ks:g.ke]
w = @view f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh]
xh = @view g.xh[g.is:g.ie+1]
z = @view g.z[g.ks:g.ke]
zh = @view g.zh[g.ks:g.keh]
rand2d = rand(g.imax, g.jmax)
rand2d .-= mean(rand2d)
s[:, :, 1] .+= rand2d[:, :]
f.s_gradbot[:, :] .= - 0.1 / settings[1]["fields"]["visc"]
f.s_gradtop[:, :] .= 0.003


@tullio s[i, j, k] += 0.003 * z[k]

xsize_max = settings[1]["grid"]["xsize"]
zsize_max = settings[1]["grid"]["zsize"]
@tullio u[i, j, k] = 10 * (0.9 + 0.1*xh[i] / xsize_max) * (z[k] / zsize_max)
# @tullio u[i, j, k] = 10 * ((xsize_max - xh[i]) / xsize_max) * (z[k] / zsize_max)

for k in 2:g.ktoth
    for i in 1:g.itot
        w[i, :, k] .= w[i, :, k-1] .- (u[i+1, :, k-1] .- u[i, :, k-1]) .* g.dz[k-1] ./ g.dx
    end
end

# w_mean = zeros(1, g.kcells)
# dudx = 1 / g.xsize .* g.z[:] ./ g.zsize
# dw = - dudx[:] .* g.dz[:]
# 
# for k in g.ks+1:g.keh
#     w_mean[1, k] = w_mean[1, k-1] + dw[k-1]
#     println("$k, $(w_mean[1, k]), $(f.w[2, 2, k])")
# end


u_gradtop = @view f.u_gradtop[g.is:g.ie+1, g.js:g.je]
@tullio u_gradtop[i, j] = (10 / zsize_max) * (xh[i] / xsize_max)


# Set the reference profile.
s_ref = @view f.s_ref[g.ks:g.ke]
s_ref[:] .= 0.003 * z[:]


## Save the restart files.
save_model(m)

