## User input.
npx = 2; npy = 2
itot = 512; jtot = 512; ktot = 512

imax = itot ÷ npx; jmax = jtot ÷ npy


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


function write_3d_hdf(a)
    fid = h5open("test_mpi.h5", "w", commxy, info)
    data_all = create_dataset(fid, "a", datatype(eltype(a)), dataspace((itot, jtot, ktot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    is = id_x*imax + 1; ie = (id_x+1)*imax
    js = id_y*jmax + 1; je = (id_y+1)*jmax
    data_all[is:ie, js:je, :] = a[:, :, :]
    MPI.Barrier(commxy)
    close(fid)
end


function write_3d_mpi(a)
    is = id_x*imax + 1; ie = (id_x+1)*imax
    js = id_y*jmax + 1; je = (id_y+1)*jmax
    subarray = MPI.Types.create_subarray((itot, jtot, ktot), (imax, jmax, ktot), (is-1, js-1, 0), MPI.Datatype(eltype(a)))
    MPI.Types.commit!(subarray)
    fid = MPI.File.open(commxy, "test_mpi.bin"; read=false, create=true, write=true, uniqueopen=true)
    MPI.File.set_view!(fid, 0, MPI.Datatype(eltype(a)), subarray, "native")
    MPI.File.write_all(fid, a[:, :, :])
    close(fid)
end


# function write_2d(a)
#     fid = h5open("test_mpi_2d.h5", "w", commx, info)
#     if id_y == 0
#         data_all = create_dataset(fid, "a2d", datatype(eltype(a)), dataspace((itot, jtot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
#         is = id_x*imax + 1; ie = (id_x+1)*imax
#         js = id_y*jmax + 1; je = (id_y+1)*jmax
#         data_all[is:ie, js:je] = a[:, :, 1]
#     end
#     MPI.Barrier(commx)
#     close(fid)
# end


for i in 1:10
    @time write_3d_hdf(a)
end

print("HDF done on id = $id\n")


# for i in 1:10
#     @time write_3d_mpi(a)
# end

# @time write_3d_mpi(a)
# @time write_3d_mpi(a)
# @time write_2d(a)
# @time write_2d(a)
# @time write_2d(a)

MPI.Finalize()
