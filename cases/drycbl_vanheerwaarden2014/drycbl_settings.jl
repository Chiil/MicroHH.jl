## Grid generation function.
function make_grid()
    # set the height (ktot = 512)
    # kmax = 512
    # dn = 1. / kmax
    
    # n = range(dn, 1-dn, length=kmax) |> collect
    
    # nloc1 = 80*dn
    # nbuf1 = 16*dn
    
    # nloc2 = 512*dn
    # nbuf2 = 72*dn
    
    # dz1 = 0.001
    # dz2 = 0.002
    # dz3 = 0.016

    kmax = 1024
    dn = 1. / kmax
    
    n = range(dn, 1-dn, length=kmax) |> collect
    
    nloc1 = 150*dn
    nbuf1 = 32*dn
    
    nloc2 = 1024*dn
    nbuf2 = 192*dn
    
    dz1 = 0.0004
    dz2 = 0.0009765625
    dz3 = 0.008
    
    dzdn1 = dz1/dn
    dzdn2 = dz2/dn
    dzdn3 = dz3/dn
    
    dzdn = similar(n)
    @. dzdn = dzdn1 + 0.5*(dzdn2-dzdn1)*(1 + tanh((n-nloc1)/nbuf1)) + 0.5*(dzdn3-dzdn2)*(1 + tanh((n-nloc2)/nbuf2))
    
    dz = dzdn*dn
    
    z       = similar(dz)
    stretch = similar(dz)
    
    z[1] = 0.5*dz[1]
    stretch[1] = 1.
    
    for k in 2:kmax
        z[k] = z[k-1] + 0.5*(dz[k-1]+dz[k])
        stretch[k] = dz[k]/dz[k-1]
    end
    
    zsize = z[kmax] + 0.5*dz[kmax]

    return z, zsize
end

z, zsize = make_grid()


## Settings.
float_type = Float32

settings_grid = Dict{String, Any}(
    "itot" => 1024,
    "jtot" => 1,
    "ktot" => 1024,

    "xsize" => 1.,
    "ysize" => 1/1024,
    "zsize" => zsize)

settings_fields = Dict(
    "visc" => 4e-5,
    "alpha" => 1)

settings_boundary = Dict(
    "mom_bot_type" => "Dirichlet",
    "mom_top_type" => "Neumann",
    "s_bot_type" => "Neumann",
    "s_top_type" => "Neumann")

settings_timeloop = Dict(
    "start_time" => 0.,
    "end_time" => 10.,
    "save_time" => 1.,
    "check_time" => 0.01,
    "dt" => 0.002)

settings_multidomain = Dict(
    "enable_nudge" => false)

settings_d01 = Dict(
    "grid" => settings_grid,
    "fields" => settings_fields,
    "boundary" => settings_boundary,
    "timeloop" => settings_timeloop,
    "multidomain" => settings_multidomain)

settings_d01["grid"]["z"] = z

settings = [ settings_d01 ]
