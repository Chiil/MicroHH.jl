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

    x::Vector{Float64}
    xh::Vector{Float64}
    y::Vector{Float64}
    yh::Vector{Float64}
    z::Vector{Float64}
    zh::Vector{Float64}

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

    z_nogc::Vector{Real} = d["z"]

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

    dxi = 1. / dx
    dyi = 1. / dy

    x = collect(range(0.5-igc, length=icells)) .* dx
    y = collect(range(0.5-jgc, length=jcells)) .* dy
    xh = collect(range(-igc, length=icells)) .* dx
    yh = collect(range(-jgc, length=jcells)) .* dy

    z = zeros(kcells)
    z[ks:ke] = z_nogc[:]
    z[ks-1] = -z[ks]
    z[ke+1] = 2*zsize - z[ke]

    zh = zeros(kcells)
    zh[ks] = 0.
    zh[keh] = zsize
    zh[ks+1:ke] .= 0.5*(z[ks:ke-1] + z[ks+1:ke])
    zh[ks-1] = -zh[ks+1]

    dzh = zeros(kcells)
    dzh[2:kcells] .= z[2:kcells] .- z[1:kcells-1]
    dzhi = 1. ./ dzh[:]
    dzh[ks-1] = dzh[ks+1]
    dzhi[ks-1] = dzhi[ks+1]

    dz = zeros(kcells)
    dz[ks:ke] .= zh[ks+1:keh] .- zh[ks:ke]
    dzi = 1. ./ dz[:]
    dz[ks-1] = dz[ks]
    dz[keh] = dz[ke]
    dzi[ks-1] = dzi[ks]
    dzi[keh] = dzi[ke]

    dzi = 1. ./ dz[:]

    g = Grid(
        itot, jtot, ktot, ktoth,
        xsize, ysize, zsize,
        igc, jgc, kgc,
        icells, jcells, kcells,
        is, js, ks, ie, je, ke, keh,
        x, xh, y, yh, z, zh,
        dx, dy, dz, dzh,
        dxi, dyi, dzi, dzhi)
end
