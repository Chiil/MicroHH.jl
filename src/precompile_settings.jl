## Grid generation function.
function make_grid(zsize, ktot)
    dz = zsize / ktot
    z = range(0.5*dz, step=dz, length=ktot) |> collect
    return z
end


## Settings.
float_type = Float32

settings_grid = Dict{String, Any}(
    "itot" => 32,
    "jtot" => 32,
    "ktot" => 32,

    "xsize" => 3200.,
    "ysize" => 3200.,
    "zsize" => 3200.)

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
    "end_time" => 3600.,
    "save_time" => 900.,
    "check_time" => 5.,
    "dt" => 5.)

settings_multidomain = Dict(
    "enable_nudge" => false)

settings_d01 = Dict(
    "grid" => settings_grid,
    "fields" => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop,
    "multidomain" => settings_multidomain)

settings_d01["grid"]["z"] = make_grid(settings_d01["grid"]["zsize"], settings_d01["grid"]["ktot"])

settings = [ settings_d01 ]
