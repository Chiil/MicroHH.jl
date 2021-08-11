using MicroHH

settings_grid = Dict(
    "itot" => 32,
    "jtot" => 32,
    "ktot" => 32,

    "xsize" => 3200.,
    "ysize" => 3200.,
    "zsize" => 3200.)

settings_fields = Dict(
    "visc" => 1.)

settings_time = Dict(
    "start_time" => 0.,
    "end_time" => 7200.,
    "dt" => 5. )

settings = Dict(
    "grid"   => settings_grid,
    "fields" => settings_fields,
    "time"   => settings_time)

grid = MicroHH.Grid(settings["grid"])
fields = MicroHH.Fields(grid, settings["fields"])

MicroHH.advection_diffusion!(fields, grid)
