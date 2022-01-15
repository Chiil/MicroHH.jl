## User input.
npx = 32; npy = 64
itot = 2048; jtot = 2048; ktot = 1024

imax = itot ÷ npx; jmax = jtot ÷ npy; kmax = ktot ÷ npx


## Init MPI and create grid.
using MPI
using LoopVectorization

MPI.Init()

dims = [npy, npx]; periodic = [1, 1]; reorder = true
commxy = MPI.Cart_create(MPI.COMM_WORLD, dims, periodic, reorder)
mpiid = MPI.Comm_rank(commxy)
commx = MPI.Cart_sub(commxy, [false, true])

data = ones(Int, imax, jmax, ktot) * mpiid
data_new = Array{Int}(undef, itot, jmax, kmax)

function transpose_zx(data_new, data)
    sendbuf = Array{Int}(undef, imax, jmax, kmax, npx)
    recvbuf = Array{Int}(undef, imax, jmax, kmax, npx)
    
    # Load the buffer.
    for i in 1:npx
        ks = (i-1)*kmax + 1
        ke = i*kmax
        @turbo sendbuf[:, :, :, i] .= data[:, :, ks:ke]
    end
    
    # Communicate data.
    MPI.Alltoall!(MPI.UBuffer(sendbuf, imax*jmax*kmax), MPI.UBuffer(recvbuf, imax*jmax*kmax), commx)
    
    # Unload the buffer.
    for i in 1:npx
        is = (i-1)*imax + 1
        ie = i*imax
        @turbo data_new[is:ie, :, :] .= recvbuf[:, :, :, i]
    end
end

transpose_zx(data_new, data)

for n in 1:20
    dt = @elapsed transpose_zx(data_new, data)
    if mpiid == 0 println("Elapsed: $dt (s)") end
end

# print("$mpiid: $(data_new[:, 1, 1])\n")

MPI.Finalize()
