## User input.
npx = 2; npy = 2
itot = 512; jtot = 512; ktot = 512

imax = itot ÷ npx; jmax = jtot ÷ npy; kblock = ktot ÷ npx


## Init MPI and create grid.
using MPI
using HDF5
using LoopVectorization

MPI.Init()

dims = [npy, npx]; periodic = [1, 1]; reorder = true
commxy = MPI.Cart_create(MPI.COMM_WORLD, dims, periodic, reorder)
commx = MPI.Cart_sub(commxy, [false, true])
commy = MPI.Cart_sub(commxy, [true, false])

info = MPI.Info()
id = MPI.Comm_rank(commxy)
id_x = MPI.Comm_rank(commx)
id_y = MPI.Comm_rank(commy)

print("$id: has parallel HDF $(HDF5.has_parallel())\n")

a = ones(imax, jmax, ktot) * id


## Transpose functions for parallel code.
function transpose_xz(data_out, data)
    if data_out === data
        return
    elseif size(data_out) == size(data)
        @tturbo data_out[:, :, :] = data[:, :, :]
    else
        sendbuf = reshape(similar(data), (imax, jmax, kblock, npx))
        recvbuf = similar(sendbuf)

        # Load the buffer.
        for i in 1:npx
            is = (i-1)*imax + 1; ie = i*imax
            @tturbo sendbuf[:, :, :, i] .= data[is:ie, :, :]
        end

        # Communicate data.
        message_size = imax*jmax*kblock
        MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), commx)

        # Unload the buffer.
        for i in 1:npx
            ks = (i-1)*kblock + 1; ke = i*kblock
            @tturbo data_out[:, :, ks:ke] .= recvbuf[:, :, :, i]
        end
    end
end


function read_3d_hdf(a)
    js = id_y*jmax + 1; je = (id_y+1)*jmax
    ks = id_x*kblock + 1; ke = (id_x+1)*kblock
    fid = h5open("test_mpi.h5", commxy, info)
    aid = read(fid, "a", dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    ax = reshape(a, (itot, jmax, kblock))
    ax[:, :, :] = aid[:, js:je, ks:ke]
    MPI.Barrier(commxy)
    transpose_xz(a, ax)
    close(fid)
end

for i in 1:10
    @time read_3d_hdf(a)
end

MPI.Finalize()
