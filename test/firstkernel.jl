using MicroHH

settings_grid = Dict(
    "itot" => 32,
    "jtot" => 32,
    "ktot" => 32,

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

settings_time = Dict(
    "start_time" => 0.,
    "end_time" => 7200.,
    "dt" => 5. )

settings = Dict(
    "grid"     => settings_grid,
    "fields"   => settings_fields,
    "time"     => settings_time,
    "dynamics" => settings_dynamics)

model = Model(settings)

