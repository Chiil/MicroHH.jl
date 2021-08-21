## Settings.
settings_grid = Dict{String, Any}(
    "itot" => 48,
    "jtot" => 48,
    "ktot" => 48,

    "xsize" => 3200.,
    "ysize" => 3200.,
    "zsize" => 3200.)

zsize = settings_grid["zsize"]
ktot = settings_grid["ktot"]
dz = zsize / ktot
z = range(0.5*dz, step=dz, length=ktot) |> collect
settings_grid["z"] = z

settings_fields = Dict(
    "visc" => 5.,
    "alpha" => 9.81 / 300)

settings_boundary = Dict(
    "mom_bot_type" => "Dirichlet",
    "mom_top_type" => "Neumann",
    "s_bot_type" => "Neumann",
    "s_top_type" => "Neumann")

settings_timeloop = Dict(
    "start_time" => 0.,
    "end_time" => 7200.,
    "save_time" => 100.,
    "check_time" => 100.,
    "dt" => 5.)

settings_d01 = Dict(
    "grid"     => settings_grid,
    "fields"   => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop)

settings = [ settings_d01 ]
