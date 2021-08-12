using MicroHH

settings_grid = Dict(
    "itot" => 64,
    "jtot" => 64,
    "ktot" => 64,

    "xsize" => 3200.,
    "ysize" => 3200.,
    "zsize" => 3200.)

settings_dynamics = Dict(
    "swadvec" => true,
    "swdiff" => true,
    "swbuoy" => true,
    "swcoriolis" => true)

settings_fields = Dict(
    "visc" => 1.)

settings_timeloop = Dict(
    "start_time" => 0.,
    "end_time" => 7200.,
    "dt" => 5. )

settings = Dict(
    "grid"     => settings_grid,
    "fields"   => settings_fields,
    "timeloop" => settings_timeloop,
    "dynamics" => settings_dynamics)

model = Model(settings)

in_progress = true
while in_progress
    global in_progress = step_model(model)
end
