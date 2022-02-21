## Grid generation function.
function make_grid(zsize, ktot)
    dz = zsize / ktot
    z = range(0.5*dz, step=dz, length=ktot) |> collect
    return z
end


## Settings.
float_type = Float64

settings_grid = Dict{String, Any}(
    "itot" => 128*32,
    "jtot" => 64*64,
    "ktot" => 1024,

    "xsize" => 6400/4096 * 128*32,
    "ysize" => 6400/4096 * 64*64,
    "zsize" => 3200)

settings_fields = Dict(
    "visc" => 1.,
    "alpha" => 9.81 / 300)

settings_boundary = Dict(
    "mom_bot_type" => "Dirichlet",
    "mom_top_type" => "Neumann",
    "s_bot_type" => "Neumann",
    "s_top_type" => "Neumann")

settings_timeloop = Dict(
    "start_time" => 0.,
    "end_time" => 25.,
    "save_time" => 900.,
    "check_time" => 5.0,
    "dt" => 0.5)

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
