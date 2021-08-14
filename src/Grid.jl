struct Grid
    # Specified by user.
    itot::Int64
    jtot::Int64
    ktot::Int64
    ktoth::Int64

    xsize::Float64
    ysize::Float64
    zsize::Float64

    igc::Int64
    jgc::Int64
    kgc::Int64

    # Calculated.
    icells::Int64
    jcells::Int64
    kcells::Int64

    is::Int64
    js::Int64
    ks::Int64

    ie::Int64
    je::Int64
    ke::Int64
    keh::Int64

    dx::Float64
    dy::Float64
    dz::Vector{Float64}
    dzh::Vector{Float64}

    dxi::Float64
    dyi::Float64
    dzi::Vector{Float64}
    dzhi::Vector{Float64}
end

function Grid(d::Dict)
    itot = d["itot"]
    jtot = d["jtot"]
    ktot = d["ktot"]
    ktoth = ktot+1

    xsize = d["xsize"]
    ysize = d["ysize"]
    zsize = d["zsize"]

    igc = 1
    jgc = 1
    kgc = 1

    icells = itot + 2*igc
    jcells = jtot + 2*jgc
    kcells = ktot + 2*kgc

    is = igc + 1
    js = jgc + 1
    ks = kgc + 1

    ie = igc + itot
    je = jgc + jtot
    ke = kgc + ktot
    keh = ke + 1

    dx = xsize / itot
    dy = ysize / jtot
    dz = zeros(kcells)
    dzh = zeros(kcells)
    dz[:] .= zsize / ktot
    dzh[:] .= zsize / ktot

    dxi = 1. / dx
    dyi = 1. / dy
    dzi = 1. ./ dz[:]
    dzhi = 1. ./ dzh[:]

    g = Grid(
        itot, jtot, ktot, ktoth,
        xsize, ysize, zsize,
        igc, jgc, kgc,
        icells, jcells, kcells,
        is, js, ks, ie, je, ke, keh,
        dx, dy, dz, dzh,
        dxi, dyi, dzi, dzhi)
end
