using MicroHH

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
    "end_time" => 7200.,
    "dt" => 5. )

settings = Dict(
    "grid"     => settings_grid,
    "fields"   => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop)

model = Model(settings)

prepare_model!(model)

in_progress = true
while in_progress
    duration = @timed global in_progress = step_model!(model)
    println(duration.time, " ", model.timeloop.time)
end

