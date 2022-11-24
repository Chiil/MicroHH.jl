struct Grid{TF <: Union{Float32, Float64}}
    # Specified by user.
    itot::Int64
    jtot::Int64
    ktot::Int64
    ktoth::Int64

    xsize::TF
    ysize::TF
    zsize::TF

    xoffset::TF
    yoffset::TF
    zoffset::TF

    igc::Int64
    jgc::Int64
    kgc::Int64

    # Calculated.
    imax::Int64
    jmax::Int64

    iblock::Int64
    jblock::Int64
    kblock::Int64

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

    x::Vector{TF}
    xh::Vector{TF}
    y::Vector{TF}
    yh::Vector{TF}
    z::Vector{TF}
    zh::Vector{TF}

    dx::TF
    dy::TF
    dz::Vector{TF}
    dzh::Vector{TF}

    dxi::TF
    dyi::TF
    dzi::Vector{TF}
    dzhi::Vector{TF}
end

function Grid(d::Dict, p::Parallel, TF)
    itot = d["itot"]
    jtot = d["jtot"]
    ktot = d["ktot"]
    ktoth = ktot+1

    xsize = d["xsize"]
    ysize = d["ysize"]
    zsize = d["zsize"]

    xoffset = d["xoffset"]
    yoffset = d["yoffset"]
    zoffset = d["zoffset"]

    z_nogc::Vector{Real} = d["z"]

    igc = 1
    jgc = 1
    kgc = 1

    imax = itot ÷ p.npx
    jmax = jtot ÷ p.npy

    iblock = itot ÷ p.npy
    jblock = jtot ÷ p.npx
    kblock = ktot ÷ p.npx

    icells = imax + 2*igc
    jcells = jmax + 2*jgc
    kcells = ktot + 2*kgc

    is = igc + 1
    js = jgc + 1
    ks = kgc + 1

    ie = igc + imax
    je = jgc + jmax
    ke = kgc + ktot
    keh = ke + 1

    dx = xsize / itot
    dy = ysize / jtot

    dxi = 1. / dx
    dyi = 1. / dy

    x = (collect(range(0.5-igc, length=icells)) .+ p.id_x*imax) .* dx .+ xoffset
    y = (collect(range(0.5-jgc, length=jcells)) .+ p.id_y*jmax) .* dy .+ yoffset
    xh = (collect(range(-igc, length=icells)) .+ p.id_x*imax) .* dx .+ xoffset
    yh = (collect(range(-jgc, length=jcells)) .+ p.id_y*jmax) .* dy .+ yoffset

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

    z[:] .+= zoffset
    zh[:] .+= zoffset

    g = Grid{TF}(
        itot, jtot, ktot, ktoth,
        xsize, ysize, zsize,
        xoffset, yoffset, zoffset,
        igc, jgc, kgc,
        imax, jmax, iblock, jblock, kblock,
        icells, jcells, kcells,
        is, js, ks, ie, je, ke, keh,
        x, xh, y, yh, z, zh,
        dx, dy, dz, dzh,
        dxi, dyi, dzi, dzhi)
end
