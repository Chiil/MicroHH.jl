## Grid generation function.
function make_grid(zsize, ktot)
    dz = zsize / ktot
    z = range(0.5*dz, step=dz, length=ktot) |> collect
    return z
end


## Settings.
settings_parallel = Dict(
    "npx" => 1,
    "npy" => 1)

settings_grid = Dict{String, Any}(
    "itot" => 256,
    "jtot" => 256,
    "ktot" => 256,

    "xsize" => 3200.,
    "ysize" => 3200.,
    "zsize" => 3200.)

settings_fields = Dict(
    "visc" => 5.,
    "alpha" => 9.81 / 300,
    "scalars" => ["smoke"])

settings_boundary = Dict(
    "mom_bot_type" => "Dirichlet",
    "mom_top_type" => "Neumann",
    "s_bot_type" => "Neumann",
    "s_top_type" => "Neumann")

settings_timeloop = Dict(
    "start_time" => 0.,
    "end_time" => 3600.,
    "save_time" => 900.,
    "check_time" => 2.,
    "dt" => 2.)

settings_d01 = Dict(
    "parallel" => settings_parallel,
    "grid" => settings_grid,
    "fields" => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop)

settings_d01["grid"]["z"] = make_grid(settings_d01["grid"]["zsize"], settings_d01["grid"]["ktot"])

settings = [ settings_d01 ]
