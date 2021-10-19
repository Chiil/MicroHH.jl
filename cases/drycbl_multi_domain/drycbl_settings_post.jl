function make_grid(zsize, ktot)
    dz = zsize / ktot
    z = range(0.5*dz, step=dz, length=ktot) |> collect
    return z
end

## Settings domain 1.
settings_grid = Dict{String, Any}(
    "itot" => 128,
    "jtot" => 128,
    "ktot" => 128,

    "xsize" => 3200.,
    "ysize" => 3200.,
    "zsize" => 3200.)

settings_fields = Dict(
    "visc" => 5.,
    "alpha" => 9.81 / 300.)

settings_boundary = Dict(
    "mom_bot_type" => "Dirichlet",
    "mom_top_type" => "Neumann",
    "s_bot_type" => "Neumann",
    "s_top_type" => "Neumann")

settings_timeloop = Dict(
    "start_time" => 7200.,
    "end_time" => 7200.,
    "save_time" => 1800.,
    "check_time" => 50.,
    "dt" => 5.)

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
settings_d02["grid"]["itot"] = 384
settings_d02["grid"]["jtot"] = 384
settings_d02["grid"]["ktot"] = 384
settings_d02["fields"]["visc"] /= 3^(4/3)
settings_d02["timeloop"]["dt"] /= 3.
settings_d02["multidomain"]["enable_nudge"] = true
settings_d02["multidomain"]["nudge_time"] = 600

settings_d02["grid"]["z"] = make_grid(settings_d02["grid"]["zsize"], settings_d02["grid"]["ktot"])


## Create the vector of Dicts, one per domain.
settings = [ settings_d01, settings_d02 ]
