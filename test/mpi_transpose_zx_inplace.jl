## User input.
npx = 4; npy = 2
itot = 256; jtot = 192; ktot = 128

imax = itot ÷ npx; jmax = jtot ÷ npy; kmax = ktot ÷ npx


## Init MPI and create grid.
using MPI
using BenchmarkTools

MPI.Init()

dims = [npy, npx]; periodic = [1, 1]; reorder = true
commxy = MPI.Cart_create(MPI.COMM_WORLD, dims, periodic, reorder)
mpiid = MPI.Comm_rank(commxy)
commx = MPI.Cart_sub(commxy, [false, true])

data = ones(Int, imax, jmax, ktot) * mpiid
buffer = Array{Int}(undef, imax, jmax, kmax, npx)

function transpose_zx(data, buffer)
    # Load the buffer.
    for i in 1:npx
        ks = (i-1)*kmax + 1
        ke = i*kmax
        buffer[:, :, :, i] .= data[:, :, ks:ke]
    end
    
    # Communicate data.
    MPI.Alltoall!(MPI.UBuffer(buffer, imax*jmax*kmax), commx)

    # Reshape the data array.
    data_xz = reshape(data, (itot, jmax, kmax))
    
    # Unload the buffer.
    for i in 1:npx
        is = (i-1)*imax + 1
        ie = i*imax
        data_xz[is:ie, :, :] .= buffer[:, :, :, i]
    end

    data_xz
end

data_xz = transpose_zx(data, buffer)

for n in 1:10
    dt = @elapsed transpose_zx(data, buffer)
    if mpiid == 0 println("Elapsed: $dt (s)") end
end

# print("$mpiid: $(data_xz[:, 1, 1])\n")

MPI.Finalize()
