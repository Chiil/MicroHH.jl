## Loading packages.
using MicroHH
using PyPlot
using MicroHH.StencilBuilder
using LoopVectorization


## Loading settings.
include("drycbl_settings.jl")
settings[1]["timeloop"]["start_time"] = 5400
settings[2]["timeloop"]["start_time"] = 5400


## Initialize the model.
n_domains = 2
m = Model("drycbl", n_domains, settings, float_type)
load_model!(m)


## Plot the output.
g1 = m.grid[1]; f1 = m.fields[1]
x1 = @view g1.xh[g1.is:g1.ie+1]
y1 = @view g1.yh[g1.js:g1.je+1]
z1 = @view g1.z[g1.ks:g1.ke]
println("z1: ", z1[8])
s1 = @view f1.s[g1.is:g1.ie, g1.js:g1.je, 8]

g2 = m.grid[2]; f2 = m.fields[2]
x2 = @view g2.xh[g2.is:g2.ie+1]
y2 = @view g2.yh[g2.js:g2.je+1]
z2 = @view g2.z[g2.ks:g2.ke]
println("z2: ", z2[23])
s2 = @view f2.s[g2.is:g2.ie, g2.js:g2.je, 23]

plt.figure()
plt.pcolormesh(x1, y1, s1', cmap=plt.cm.turbo, vmin=1.75, vmax=2.5)
plt.pcolormesh(x2, y2, s2', cmap=plt.cm.turbo, vmin=1.75, vmax=2.5)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.gca().set_aspect("equal")
plt.tight_layout()
display(plt.gcf())
plt.show()
