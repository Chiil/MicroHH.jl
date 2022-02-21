## User input.
npx = 32; npy = 64
itot = 4096; jtot = 4096; ktot = 1024

imax = itot ÷ npx; jmax = jtot ÷ npy; kmax = ktot ÷ npx


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

a = ones(imax, jmax, ktot) * id

function write_3d_hdf(a)
    is = id_x*imax + 1; ie = (id_x+1)*imax
    js = id_y*jmax + 1; je = (id_y+1)*jmax

    fid = h5open("test_mpi.h5", "w", commxy, info)

    # aid = create_dataset(fid, "a", datatype(a), dataspace(itot, jtot, ktot), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    aid = create_dataset(fid, "a", datatype(a), dataspace(itot, jtot, ktot), chunk=(imax, jmax, ktot), dxpl_mpio=HDF5.H5FD_MPIO_COLLECTIVE)
    aid[is:ie, js:je, :] = a[:, :, :]
    close(aid)

    close(fid)
end

function read_3d_hdf(a)
    is = id_x*imax + 1; ie = (id_x+1)*imax
    js = id_y*jmax + 1; je = (id_y+1)*jmax

    fid = h5open("test_mpi.h5", commxy, info)

    dapl = create_property(HDF5.H5P_DATASET_ACCESS)
    dxpl = create_property(HDF5.H5P_DATASET_XFER)
    dxpl[:dxpl_mpio] = HDF5.H5FD_MPIO_COLLECTIVE
    aid = open_dataset(fid, "a", dapl, dxpl)
    a[:, :, :] = aid[is:ie, js:je, :]
    close(aid)

    close(fid)
end

for i in 1:3
    if id == 0
        @time write_3d_hdf(a)
    else
        write_3d_hdf(a)
    end
end

for i in 1:3
    if id == 0
        @time read_3d_hdf(a)
    else
        read_3d_hdf(a)
    end
end

MPI.Finalize()
