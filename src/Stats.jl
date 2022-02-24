module Stats

export open_stats, add_dimension!, add_field!, add_record!, add_time_record!
using HDF5


function open_stats(filename::String)
    fid = h5open(filename, "cw")
    return fid
end


function add_dimension!(fid, d_name, d)
    if ndims(d) != 1
        @error "Dimension $d_name is not 1D"
    end

    if haskey(fid, d_name)
        @error "Dimension $d_name already exists"
    else
        write(fid, d_name, d[:])
        HDF5.h5ds_set_scale(fid[d_name], d_name)
    end
end


function add_field!(fid, a, a_name, dims)
    if haskey(fid, a_name)
        @error "Field already exists"
    else
        # Create variable.
        @info "Saving variable $a_name"
        aid = create_dataset(fid, a_name, datatype(eltype(a)), dataspace(size(a)...))

        # Add dimensions to variable.
        for i in 1:length(dims)
            index_c = length(dims)-i
            HDF5.h5ds_attach_scale(fid[a_name], fid[dims[i]], index_c)
        end

        # Save the first data.
        r = [ 1:lastindex(a, i) for i in 1:ndims(a) ]
        aid[r...] = a[r...]

        close(aid)
    end
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

        if length(dims) == 0
            HDF5.h5ds_set_scale(fid["time"], "time")
        else
            # Add dimensions to variable.
            # Append time as it is a record.
            dims = (dims..., "time")
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
        @error "Field already exist"
    end

    close(aid)
    close(iterid)
end


add_time_record!(fid, time, iter) = add_record!(fid, time, "time", (), iter)


end
