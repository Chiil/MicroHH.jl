## Grid generation function.
function make_grid(zsize, ktot)
    dz = zsize / ktot
    z = range(0.5*dz, step=dz, length=ktot) |> collect
    return z
end


## Settings.
float_type = Float64

settings_grid = Dict{String, Any}(
    "itot" => 128,
    "jtot" => 1,
    "ktot" => 128,

    "xsize" => 1.,
    "ysize" => 1. / 256.,
    "zsize" => 1.,

    "xoffset" => 0.,
    "yoffset" => 0.,
    "zoffset" => 0.)

settings_fields = Dict(
    "visc" => 1. / 100.,
    "alpha" => 0.)

settings_boundary = Dict(
    "mom_bot_type" => "Dirichlet",
    "mom_top_type" => "Dirichlet",
    "s_bot_type" => "Neumann",
    "s_top_type" => "Neumann")

settings_timeloop = Dict(
    "start_time" => 0.,
    "end_time" => 10.,
    "save_time" => 1.,
    "check_time" => 0.1,
    "dt" => 0.001)

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