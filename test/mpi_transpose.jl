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
data_new = Array{Int}(undef, itot, jmax, kmax)

function transpose_zx(data_new, data)
    sendbuf = Array{Int}(undef, imax, jmax, kmax, npx)
    recvbuf = Array{Int}(undef, imax, jmax, kmax, npx)
    
    # Load the buffer
    for i in 1:npx
        ks = (i-1)*kmax + 1
        ke = i*kmax
        sendbuf[:, :, :, i] .= data[:, :, ks:ke]
    end
    
    # Transpose zx
    # reqs = Vector{MPI.Request}(undef, 2*npx)
    # for i in 1:npx
    #     tag = 1
    #     r = @view recvbuf[:, :, :, i]; s = @view sendbuf[:, :, :, i]
    #     reqs[2*(i-1)+1] = MPI.Irecv!(r, i-1, tag, commx)
    #     reqs[2*(i-1)+2] = MPI.Isend( s, i-1, tag, commx)
    # end
    # MPI.Waitall!(reqs)
    MPI.Alltoall!(MPI.UBuffer(sendbuf, imax*jmax*kmax), MPI.UBuffer(recvbuf, imax*jmax*kmax), commx)
    
    # Unload the buffer
    for i in 1:npx
        is = (i-1)*imax + 1
        ie = i*imax
        data_new[is:ie, :, :] .= recvbuf[:, :, :, i]
    end
end

transpose_zx(data_new, data)

for n in 1:10
    dt = @elapsed transpose_zx(data_new, data)
    if mpiid == 0 println("Elapsed: $dt (s)") end
end

print("$mpiid: $(data_new[:, 1, 1])\n")

MPI.Finalize()
