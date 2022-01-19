## User input.
npx = 2; npy = 4
itot = 4; jtot = 8; ktot = 6

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
fid = h5open("test_mpi.h5", "w", commxy, info)
data_all = create_dataset(fid, "a", datatype(eltype(a)), dataspace((itot, jtot, ktot)), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
is = id_x*imax + 1; ie = (id_x+1)*imax
js = id_y*jmax + 1; je = (id_y+1)*jmax
print("$id, $id_x, $id_y, $is, $ie, $js, $je\n")
data_all[is:ie, js:je, :] = a[:, :, :]

close(fid)
MPI.Finalize()
