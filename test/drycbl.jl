## Loading packages.
using MicroHH
using Tullio
using HDF5
using Statistics


## Model settings.
settings_grid = Dict(
    "itot" => 256,
    "jtot" => 256,
    "ktot" => 256,

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
    "end_time" => 1800.,
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
z = range(g.dz[1]/2, g.zsize, step=g.dz[1]) |> collect
rand2d = rand(g.itot, g.jtot)
rand2d .-= mean(rand2d)
s[:, :, 1] .+= rand2d[:, :]
@tullio s[i, j, k] += 0.003 * z[k]


## Run the model.
prepare_model!(model)

in_progress = true
while in_progress
    duration = @timed global in_progress = step_model!(model)
    println(duration.time, " ", model.timeloop.time)
end


## Save some fields.
h5open("drycbl.end.h5", "w") do fid
    write(fid, "u", model.fields.u[:, :, :])
    write(fid, "v", model.fields.v[:, :, :])
    write(fid, "w", model.fields.w[:, :, :])
    write(fid, "s", model.fields.s[:, :, :])
end
