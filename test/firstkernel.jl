using MicroHH

settings = Dict(
    "itot" => 32,
    "jtot" => 32,
    "ktot" => 32,

    "xsize" => 3200.,
    "ysize" => 3200.,
    "zsize" => 3200.)

grid = MicroHH.Grid(settings)
fields = MicroHH.Fields(grid)

MicroHH.advection_diffusion!(fields, grid)
