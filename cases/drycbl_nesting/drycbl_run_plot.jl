## Loading packages.
using MicroHH
using Tullio
using Statistics
using GLMakie


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings, float_type)
load_model!(m)
in_progress = prepare_model!(m)


## Make the figure.
fig = Figure(resolution=(800, 500))

f = m.fields[1]; g = m.grid[1]
x = @view g.x[g.is:g.ie]
y = @view g.y[g.js:g.je]
z = @view g.z[g.ks:g.ke]
s_bot = @view f.s_bot[g.is:g.ie, g.js:g.je]
s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]
s_ref = reshape(f.s_ref[g.ks:g.ke], (1, 1, length(z)))

node1 = Observable(s_bot .- mean(s_bot))
node2 = Observable(s[:, 1, :] .- s_ref[:, 1, :])
h1 = heatmap(fig[1, 1], x, y, node1, colorrange=(-2, 2))
h2 = heatmap(fig[1, 2], x, z, node2, colorrange=(0, 2))
h1ref = heatmap(fig[2, 1], x, y, node1, colorrange=(-2, 2))
h2ref = heatmap(fig[2, 2], x, z, node2, colorrange=(0, 2))

if n_domains > 1
    f2 = m.fields[2]; g2 = m.grid[2]
    f2 = m.fields[2]; g2 = m.grid[2]
    x2 = @view g2.x[g2.is:g2.ie]
    y2 = @view g2.y[g2.js:g2.je]
    z2 = @view g2.z[g2.ks:g2.ke]
    s_bot2 = @view f2.s_bot[g2.is:g2.ie, g2.js:g2.je]
    s2 = @view f2.s[g2.is:g2.ie, g2.js:g2.je, g2.ks:g2.ke]
    s_ref2 = reshape(f2.s_ref[g2.ks:g2.ke], (1, 1, length(z2)))

    node3 = Observable(s_bot2 .- mean(s_bot))
    node4 = Observable(s2[:, 1, :] .- s_ref2[:, 1, :])
    h3 = heatmap!(fig[1, 1], x2, y2, node3, colorrange=(-2, 2))
    h4 = heatmap!(fig[1, 2], x2, z2, node4, colorrange=(0, 2))
end

display(fig)


## Run the model.
while in_progress
    global in_progress = step_model!(m)
    node1[] = s_bot .- mean(s_bot)
    node2[] = s[:, 1, :] .- s_ref[:, 1, :]

    if n_domains > 1
        node3[] = s_bot2 .- mean(s_bot)
        node4[] = s2[:, 1, :] .- s_ref2[:, 1, :]
    end
end
