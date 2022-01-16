## User input.
# npx = 2; npy = 4
# itot = 8; jtot = 8; ktot = 8
npx = 32; npy = 64
itot = 2048; jtot = 2048; ktot = 1024


## Init MPI and create grid.
using MPI
using LoopVectorization

imax = itot ÷ npx; jmax = jtot ÷ npy
iblock = itot ÷ npy; jblock = jtot ÷ npx; kblock = ktot ÷ npx

MPI.Init()

dims = [npy, npx]; periodic = [1, 1]; reorder = true
commxy = MPI.Cart_create(MPI.COMM_WORLD, dims, periodic, reorder)
mpiid = MPI.Comm_rank(commxy)
commx = MPI.Cart_sub(commxy, [false, true])
commy = MPI.Cart_sub(commxy, [true, false])

data_k = ones(Int, imax, jmax, ktot) * mpiid
data_i = Array{Int}(undef, itot, jmax, kblock)
data_j = Array{Int}(undef, iblock, jtot, kblock)
data_kt = Array{Int}(undef, iblock, jblock, ktot)


## Create transpose functions.
function transpose_zx(data_new, data)
    sendbuf = Array{Int}(undef, imax, jmax, kblock, npx)
    recvbuf = Array{Int}(undef, imax, jmax, kblock, npx)

    # Load the buffer.
    for i in 1:npx
        ks = (i-1)*kblock + 1; ke = i*kblock
        @tturbo sendbuf[:, :, :, i] .= data[:, :, ks:ke]
    end

    # Communicate data.
    MPI.Alltoall!(MPI.UBuffer(sendbuf, imax*jmax*kblock), MPI.UBuffer(recvbuf, imax*jmax*kblock), commx)

    # Unload the buffer.
    for i in 1:npx
        is = (i-1)*imax + 1; ie = i*imax
        @tturbo data_new[is:ie, :, :] .= recvbuf[:, :, :, i]
    end
end


function transpose_xy(data_new, data)
    sendbuf = Array{Int}(undef, iblock, jmax, kblock, npy)
    recvbuf = Array{Int}(undef, iblock, jmax, kblock, npy)

    # Load the buffer.
    for i in 1:npy
        is = (i-1)*iblock + 1; ie = i*iblock
        @tturbo sendbuf[:, :, :, i] .= data[is:ie, :, :]
    end

    # Communicate data.
    MPI.Alltoall!(MPI.UBuffer(sendbuf, iblock*jmax*kblock), MPI.UBuffer(recvbuf, iblock*jmax*kblock), commy)

    # Unload the buffer.
    for i in 1:npy
        js = (i-1)*jmax + 1; je = i*jmax
        @tturbo data_new[:, js:je, :] .= recvbuf[:, :, :, i]
    end
end


function transpose_yzt(data_new, data)
    sendbuf = Array{Int}(undef, iblock, jblock, kblock, npx)
    recvbuf = Array{Int}(undef, iblock, jblock, kblock, npx)

    # Load the buffer.
    for i in 1:npx
        js = (i-1)*jblock + 1; je = i*jblock
        @tturbo sendbuf[:, :, :, i] .= data[:, js:je, :]
    end

    # Communicate data.
    MPI.Alltoall!(MPI.UBuffer(sendbuf, iblock*jblock*kblock), MPI.UBuffer(recvbuf, iblock*jblock*kblock), commx)

    # Unload the buffer.
    for i in 1:npx
        ks = (i-1)*kblock + 1; ke = i*kblock
        @tturbo data_new[:, :, ks:ke] .= recvbuf[:, :, :, i]
    end
end


function transpose_zty(data_new, data)
    sendbuf = Array{Int}(undef, iblock, jblock, kblock, npx)
    recvbuf = Array{Int}(undef, iblock, jblock, kblock, npx)

    # Load the buffer.
    for i in 1:npx
        ks = (i-1)*kblock + 1; ke = i*kblock
        @tturbo sendbuf[:, :, :, i] .= data[:, :, ks:ke]
    end

    # Communicate data.
    MPI.Alltoall!(MPI.UBuffer(sendbuf, iblock*jblock*kblock), MPI.UBuffer(recvbuf, iblock*jblock*kblock), commx)

    # Unload the buffer.
    for i in 1:npx
        js = (i-1)*jblock + 1; je = i*jblock
        @tturbo data_new[:, js:je, :] = recvbuf[:, :, :, i]
    end
end


function transpose_yx(data_new, data)
    sendbuf = Array{Int}(undef, iblock, jmax, kblock, npy)
    recvbuf = Array{Int}(undef, iblock, jmax, kblock, npy)

    # Load the buffer.
    for i in 1:npy
        js = (i-1)*jmax + 1; je = i*jmax
        @tturbo sendbuf[:, :, :, i] .= data[:, js:je, :]
    end

    # Communicate data.
    MPI.Alltoall!(MPI.UBuffer(sendbuf, iblock*jmax*kblock), MPI.UBuffer(recvbuf, iblock*jmax*kblock), commy)

    # Unload the buffer.
    for i in 1:npy
        is = (i-1)*iblock + 1; ie = i*iblock
        @tturbo data_new[is:ie, :, :] .= recvbuf[:, :, :, i]
    end
end


function transpose_xz(data_new, data)
    sendbuf = Array{Int}(undef, imax, jmax, kblock, npx)
    recvbuf = Array{Int}(undef, imax, jmax, kblock, npx)

    # Load the buffer.
    for i in 1:npx
        is = (i-1)*imax + 1; ie = i*imax
        @tturbo sendbuf[:, :, :, i] .= data[is:ie, :, :]
    end

    # Communicate data.
    MPI.Alltoall!(MPI.UBuffer(sendbuf, imax*jmax*kblock), MPI.UBuffer(recvbuf, imax*jmax*kblock), commx)

    # Unload the buffer.
    for i in 1:npx
        ks = (i-1)*kblock + 1; ke = i*kblock
        @tturbo data_new[:, :, ks:ke] .= recvbuf[:, :, :, i]
    end
end


## Test the transposes.
if mpiid == 0 print("Testing transposes with npx = $npx, npy = $npy\n") end
MPI.Barrier(commxy)

# print("(s , k) $mpiid: $(data_k[1, 1, :])\n")
# transpose_zx(data_i, data_k)
# print("(zx, i) $mpiid: $(data_i[:, 1, 1])\n")
# transpose_xy(data_j, data_i)
# print("(xy, j) $mpiid: $(data_j[1, :, 1])\n")
# transpose_yzt(data_kt, data_j)
# print("(yz, k) $mpiid: $(data_kt[1, 1, :])\n")
# transpose_zty(data_j, data_kt)
# print("(zy, j) $mpiid: $(data_j[1, :, 1])\n")
# transpose_yx(data_i, data_j)
# print("(xy, i) $mpiid: $(data_i[:, 1, 1])\n")
# transpose_xz(data_k, data_i)
# print("(xz, k) $mpiid: $(data_k[1, 1, :])\n")


## Run a timed benchmark
transpose_zx(data_i, data_k)
transpose_xy(data_j, data_i)
transpose_yzt(data_kt, data_j)
transpose_zty(data_j, data_kt)
transpose_yx(data_i, data_j)
transpose_xz(data_k, data_i)

for n in 1:20
    MPI.Barrier(commxy)
    dt = @elapsed begin
        transpose_zx(data_i, data_k)
        transpose_xy(data_j, data_i)
        transpose_yzt(data_kt, data_j)
        transpose_zty(data_j, data_kt)
        transpose_yx(data_i, data_j)
        transpose_xz(data_k, data_i)
    end
    if mpiid == 0 println("Elapsed: $dt (s)") end
end

MPI.Barrier(commxy)
if mpiid == 0 print("Test completed.\n") end


## Close the MPI.
MPI.Finalize()
