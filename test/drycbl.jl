## Loading packages.
using MicroHH
using Tullio
using HDF5


## Model settings.
settings_grid = Dict(
    "itot" => 32,
    "jtot" => 32,
    "ktot" => 32,

    "xsize" => 3200.,
    "ysize" => 3200.,
    "zsize" => 3200.)

settings_fields = Dict(
    "visc" => 1.)

settings_boundary = Dict(
    "mom_bot_type" => "Dirichlet",
    "mom_top_type" => "Neumann",
    "s_bot_type" => "Dirichlet",
    "s_top_type" => "Neumann")

settings_timeloop = Dict(
    "start_time" => 0.,
    "end_time" => 7200.,
    "dt" => 5. )

settings = Dict(
    "grid"     => settings_grid,
    "fields"   => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop)

model = Model(settings)


## Fields initialization.
f = model.fields; g = model.grid
s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke]
u = @view f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke]
v = @view f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke]
w = @view f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh]

z = range(g.dz[1]/2, g.zsize, step=g.dz[1]) |> collect
s_rand = rand(g.itot, g.jtot) .- 0.5
s[:, :, 1] .+= s_rand[:, :]
@tullio s[i, j, k] += 0.003 * z[k]
h5open("drycbl.00000000.h5", "w") do fid
    write(fid, "u", u[:, :, :])
    write(fid, "v", v[:, :, :])
    write(fid, "w", w[:, :, :])
    write(fid, "s", s[:, :, :])
end


## Run the model.
prepare_model!(model)

in_progress = true
while in_progress
    duration = @timed global in_progress = step_model!(model)
    println(duration.time, " ", model.timeloop.time)
end

