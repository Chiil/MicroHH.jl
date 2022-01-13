using HDF5

file_id = h5open("test_unlimited.h5", "cw")

a = rand(2, 3, 4)

a_id = create_dataset(file_id, "a", Float64, ((size(a)..., 1), (size(a)..., -1)), chunk=(size(a)..., 1))
a_id[:, :, :, 1] .= a[:, :, :]

for i in 1:10
    a[:, :, :] .+= 1
    HDF5.set_extent_dims(a_id, (size(a)..., i))
    a_id[:, :, :, i] = a[:, :, :]
end

close(file_id)
