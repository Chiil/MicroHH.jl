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
        iterid = create_dataset(fid, "iter", datatype(typeof(iter)), ((1,), (-1,)), chunk=(1,))
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
end


## Save a field
time = 0.

fid = open_stats("test.h5")
add_record!(fid, time, "time", :t, 0)
add_record!(fid, time, "time", :t, 10)
add_record!(fid, time, "time", :t, 300)
close(fid)
