## Load packages.
using HDF5


## Stats functions.
function open_stats(filename::String)
    fid = h5open(filename, "cw")
    return fid
end


function add_record!(fid, a, a_name, dims, iter)
    # Check if record variables are available.
    if !haskey(fid, "iter")
        @info "Iter does not exist, incrementing size and adding iter"
        iterid = create_dataset(fid, "iter", datatype(eltype(iter)), ((1,), (-1,)), chunk=(1,))
        iterid[1] = iter
    else
        iterid = fid["iter"]

        # Check if iter exists in data
        if iter in iterid[:]
            @info "Iter already exists"
        else
            @info "Iter does not exist, incrementing size and adding iter"
            iterdims = [ size(iterid)... ]
            iterdims[end] += 1
            HDF5.set_extent_dims(iterid, tuple(iterdims...))
            iterid[end] = iter
        end
    end

    # Iter is available, continue to store data.
    if !haskey(fid, a_name)
        # Create variable.
        @info "Variable $a_name does not exist, incrementing size and adding variable"
        aid = create_dataset(fid, a_name, datatype(eltype(a)), ((size(a)..., 1,), (size(a)..., -1,)), chunk=(size(a)..., 1,))

        if dims == nothing
            HDF5.h5ds_set_scale(fid["time"], "time")
        else
            # Add dimensions to variable.
            for i in 1:length(dims)
                index_c = length(dims)-i
                HDF5.h5ds_attach_scale(fid[a_name], fid[dims[i]], index_c)
            end
        end

        # Save the first data.
        @info "Saving data..."
        r = [ 1:lastindex(a, i) for i in 1:ndims(a) ]
        aid[r..., end] = a[r...]
        return nothing
    else
        aid = fid[a_name]
    end

    if size(aid)[end] == size(iterid)[end]-1
        @info "Saving data..."
        adims = [ size(aid)... ]
        adims[end] += 1
        HDF5.set_extent_dims(aid, tuple(adims...))
        r = [ 1:lastindex(a, i) for i in 1:ndims(a) ]
        aid[r..., end] = a[r...]
    else
        @info "OHOH DATA IS ALREADY THERE!"
    end

    close(aid)
    close(iterid)
end


## Save a field
fid = open_stats("test.h5")

time = 0.
a = rand()

add_record!(fid, time, "time", nothing, 0)
add_record!(fid, a, "a", ("time",), 0)

time += 10.
a = rand()
add_record!(fid, time, "time", nothing, 10)
add_record!(fid, a, "a", ("time",), 10)

time += 300.
a = rand()
add_record!(fid, time, "time", nothing, 300)
add_record!(fid, a, "a", ("time",), 300)

close(fid)
