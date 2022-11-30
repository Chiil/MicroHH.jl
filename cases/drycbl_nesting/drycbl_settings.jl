## Grid generation function.
function make_grid(zsize, ktot)
    dz = zsize / ktot
    z = range(0.5*dz, step=dz, length=ktot) |> collect
    return z
end


## Default float type.
float_type = Float64


## Settings domain 1.
settings_grid = Dict{String, Any}(
    "itot" => 512,
    "jtot" => 1,
    "ktot" => 256,

    "xsize" => 6400.,
    "ysize" => 3200. / 256,
    "zsize" => 3200.,

    "xoffset" => 0.,
    "yoffset" => 0.,
    "zoffset" => 0.)

settings_fields = Dict(
    "visc" => 1.5,
    "alpha" => 9.81 / 300.)

settings_boundary = Dict(
    "mom_bot_type" => "Dirichlet",
    "mom_top_type" => "Neumann",
    "s_bot_type" => "Neumann",
    "s_top_type" => "Neumann")

settings_timeloop = Dict(
    "start_time" => 3600.,
    "end_time" => 7200.,
    "save_time" => 300.,
    "check_time" => 50.,
    "dt" => 1.0)

settings_multidomain = Dict{String, Any}(
    "enable_nudge" => false)

settings_d01 = Dict(
    "grid" => settings_grid,
    "fields" => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop,
    "multidomain" => settings_multidomain)

settings_d01["grid"]["z"] = make_grid(settings_d01["grid"]["zsize"], settings_d01["grid"]["ktot"])


# Settings domain 2.
settings_d02 = deepcopy(settings_d01)
settings_d02["grid"]["itot"] = 256
#settings_d02["grid"]["jtot"] = 256
settings_d02["grid"]["ktot"] = 256
settings_d02["grid"]["xsize"] = 3200.
#settings_d02["grid"]["ysize"] = 3200.
settings_d02["grid"]["zsize"] = 3200.
settings_d02["grid"]["xoffset"] = 1600.
settings_d02["grid"]["yoffset"] = 0.
settings_d02["grid"]["zoffset"] = 0.

settings_d02["grid"]["z"] = make_grid(settings_d02["grid"]["zsize"], settings_d02["grid"]["ktot"])


## Create the vector of Dicts, one per domain.
settings = [ settings_d01, settings_d02 ]
