using HDF5

filename = "test.h5"

x = 0.:10:100 |> collect
y = 0.:2:10 |> collect
z = [ i^2 + j for i in x, j in y ]

h5open(filename, "w") do fid
    fid["x"] = x[:]
    fid["y"] = y[:]
    fid["z"] = z[:, :]
    HDF5.h5ds_set_scale(fid["x"], "x")
    HDF5.h5ds_set_scale(fid["y"], "y")
    HDF5.h5ds_attach_scale(fid["z"], fid["x"], 1)
    HDF5.h5ds_attach_scale(fid["z"], fid["y"], 0)
end
