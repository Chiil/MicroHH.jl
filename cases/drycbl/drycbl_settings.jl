## Grid generation function.
function make_grid(zsize, ktot)
    dz = zsize / ktot
    z = range(0.5*dz, step=dz, length=ktot) |> collect
    return z
end


## Settings.
settings_grid = Dict{String, Any}(
    "itot" => 192,
    "jtot" => 192,
    "ktot" => 192,

    "xsize" => 3200.,
    "ysize" => 3200.,
    "zsize" => 3200.)

settings_fields = Dict(
    "visc" => 2.,
    "alpha" => 9.81 / 300)

settings_boundary = Dict(
    "mom_bot_type" => "Dirichlet",
    "mom_top_type" => "Neumann",
    "s_bot_type" => "Neumann",
    "s_top_type" => "Neumann")

settings_timeloop = Dict(
    "start_time" => 0.,
    "end_time" => 7200.,
    "save_time" => 7200.,
    "check_time" => 100.,
    "dt" => 5.)

settings_d01 = Dict(
    "grid"     => settings_grid,
    "fields"   => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop)

settings_d01["grid"]["z"] = make_grid(settings_d01["grid"]["zsize"], settings_d01["grid"]["ktot"])

settings = [ settings_d01 ]
