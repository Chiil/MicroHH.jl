function make_grid(zsize, ktot)
    dz = zsize / ktot
    z = range(0.5*dz, step=dz, length=ktot) |> collect
    return z
end

## Settings domain 1.
settings_grid = Dict{String, Any}(
    "itot" => 48,
    "jtot" => 48,
    "ktot" => 48,

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
    "start_time" => 0.,
    "end_time" => 7200.,
    "save_time" => 100000.,
    "check_time" => 100.,
    "dt" => 20.)

settings_d01 = Dict(
    "grid"     => settings_grid,
    "fields"   => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop)

settings_d01["grid"]["z"] = make_grid(settings_d01["grid"]["zsize"], settings_d01["grid"]["ktot"])


# Settings domain 2.
settings_d02 = deepcopy(settings_d01)
settings_d02["grid"]["itot"] = 96
settings_d02["grid"]["jtot"] = 96
settings_d02["grid"]["ktot"] = 96
settings_d02["timeloop"]["dt"] = 10.

settings_d02["grid"]["z"] = make_grid(settings_d02["grid"]["zsize"], settings_d02["grid"]["ktot"])


## Create the vector of Dicts, one per domain.
settings = [ settings_d01, settings_d02 ]
