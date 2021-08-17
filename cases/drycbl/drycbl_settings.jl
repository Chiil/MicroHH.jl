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
    "save_time" => 1800.,
    "dt" => 4.)

settings = Dict(
    "grid"     => settings_grid,
    "fields"   => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop)

