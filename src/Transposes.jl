## Transpose functions for serial code.
function transpose_zx!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelSerial)
    ## Check whether data and data_out point to the same memory
    if data_out === data
        return
    else
        @tturbo data_out[:, :, :] .= data
    end
end


function transpose_xy!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelSerial)
    ## Check whether data and data_out point to the same memory
    if data_out === data
        return
    else
        @tturbo data_out[:, :, :] .= data
    end
end


function transpose_yzt!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelSerial)
    ## Check whether data and data_out point to the same memory
    if data_out === data
        return
    else
        @tturbo data_out[:, :, :] .= data
    end
end


function transpose_zty!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelSerial)
    ## Check whether data and data_out point to the same memory
    if data_out === data
        return
    else
        @tturbo data_out[:, :, :] .= data
    end
end


function transpose_yx!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelSerial)
    ## Check whether data and data_out point to the same memory
    if data_out === data
        return
    else
        @tturbo data_out[:, :, :] .= data
    end
end


function transpose_xz!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelSerial)
    ## Check whether data and data_out point to the same memory
    if data_out === data
        return
    else
        @tturbo data_out[:, :, :] .= data
    end
end


## Transpose functions for parallel code.
function transpose_zx!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelDistributed)
    ## Check whether data and data_out point to the same memory
    if data_out === data
        return
    elseif size(data_out) == size(data)
        @tturbo data_out[:, :, :] .= data
    else
        sendbuf = reshape(sendbuf, (g.imax, g.jmax, g.kblock, p.npx))
        recvbuf = reshape(recvbuf, (g.imax, g.jmax, g.kblock, p.npx))

        # Load the buffer.
        for i in 1:p.npx
            ks = (i-1)*g.kblock + 1; ke = i*g.kblock
            data_view = @view data[:, :, ks:ke]
            @tturbo sendbuf[:, :, :, i] .= data_view
        end

        # Communicate data using collective
        # message_size = g.imax*g.jmax*g.kblock
        # MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commx)

        # Communicate data using point-to-point.
        reqs = Vector{MPI.Request}(undef, 2*p.npx)

        for i in 1:p.npx
            recvbuf_view = @view recvbuf[:, :, :, i];
            reqs[i] = MPI.Irecv!(recvbuf_view, i-1, 1, p.commx)
        end

        for i in 1:p.npx
            sendbuf_view = @view sendbuf[:, :, :, i]
            reqs[i+p.npx] = MPI.Isend(sendbuf_view, i-1, 1, p.commx)
        end

        MPI.Waitall!(reqs)

        # Unload the buffer.
        for i in 1:p.npx
            is = (i-1)*g.imax + 1; ie = i*g.imax
            recvbuf_view = @view recvbuf[:, :, :, i]
            @tturbo data_out[is:ie, :, :] .= recvbuf_view
        end
    end
end


function transpose_xy!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelDistributed)
    if data_out === data
        return
    elseif size(data_out) == size(data)
        @tturbo data_out[:, :, :] .= data
    else
        sendbuf = reshape(sendbuf, (g.iblock, g.jmax, g.kblock, p.npy))
        recvbuf = reshape(recvbuf, (g.iblock, g.jmax, g.kblock, p.npy))

        # Load the buffer.
        for i in 1:p.npy
            is = (i-1)*g.iblock + 1; ie = i*g.iblock
            data_view = @view data[is:ie, :, :]
            @tturbo sendbuf[:, :, :, i] .= data_view
        end

        # Communicate data.
        # message_size = g.iblock*g.jmax*g.kblock
        # MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commy)

        # Communicate data using point-to-point.
        reqs = Vector{MPI.Request}(undef, 2*p.npy)

        for i in 1:p.npy
            recvbuf_view = @view recvbuf[:, :, :, i];
            reqs[i] = MPI.Irecv!(recvbuf_view, i-1, 1, p.commy)
        end

        for i in 1:p.npy
            sendbuf_view = @view sendbuf[:, :, :, i]
            reqs[i+p.npy] = MPI.Isend(sendbuf_view, i-1, 1, p.commy)
        end

        MPI.Waitall!(reqs)

        # Unload the buffer.
        for i in 1:p.npy
            js = (i-1)*g.jmax + 1; je = i*g.jmax
            recvbuf_view = @view recvbuf[:, :, :, i]
            @tturbo data_out[:, js:je, :] .= recvbuf_view
        end
    end
end


function transpose_yzt!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelDistributed)
    if data_out === data
        return
    elseif size(data_out) == size(data)
        @tturbo data_out[:, :, :] .= data
    else
        sendbuf = reshape(sendbuf, (g.iblock, g.jblock, g.kblock, p.npx))
        recvbuf = reshape(recvbuf, (g.iblock, g.jblock, g.kblock, p.npx))

        # Load the buffer.
        for i in 1:p.npx
            js = (i-1)*g.jblock + 1; je = i*g.jblock
            data_view = @view data[:, js:je, :]
            @tturbo sendbuf[:, :, :, i] .= data_view
        end

        # Communicate data.
        # message_size = g.iblock*g.jblock*g.kblock
        # MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commx)

        # Communicate data using point-to-point.
        reqs = Vector{MPI.Request}(undef, 2*p.npx)

        for i in 1:p.npx
            recvbuf_view = @view recvbuf[:, :, :, i];
            reqs[i] = MPI.Irecv!(recvbuf_view, i-1, 1, p.commx)
        end

        for i in 1:p.npx
            sendbuf_view = @view sendbuf[:, :, :, i]
            reqs[i+p.npx] = MPI.Isend(sendbuf_view, i-1, 1, p.commx)
        end

        MPI.Waitall!(reqs)

        # Unload the buffer.
        for i in 1:p.npx
            ks = (i-1)*g.kblock + 1; ke = i*g.kblock
            recvbuf_view = @view recvbuf[:, :, :, i]
            @tturbo data_out[:, :, ks:ke] .= recvbuf_view
        end
    end
end


function transpose_zty!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelDistributed)
    if data_out === data
        return
    elseif size(data_out) == size(data)
        @tturbo data_out[:, :, :] .= data
    else
        sendbuf = reshape(sendbuf, (g.iblock, g.jblock, g.kblock, p.npx))
        recvbuf = reshape(recvbuf, (g.iblock, g.jblock, g.kblock, p.npx))

        # Load the buffer.
        for i in 1:p.npx
            ks = (i-1)*g.kblock + 1; ke = i*g.kblock
            data_view = @view data[:, :, ks:ke]
            @tturbo sendbuf[:, :, :, i] .= data_view
        end

        # Communicate data.
        # message_size = g.iblock*g.jblock*g.kblock
        # MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commx)

        # Communicate data using point-to-point.
        reqs = Vector{MPI.Request}(undef, 2*p.npx)

        for i in 1:p.npx
            recvbuf_view = @view recvbuf[:, :, :, i];
            reqs[i] = MPI.Irecv!(recvbuf_view, i-1, 1, p.commx)
        end

        for i in 1:p.npx
            sendbuf_view = @view sendbuf[:, :, :, i]
            reqs[i+p.npx] = MPI.Isend(sendbuf_view, i-1, 1, p.commx)
        end

        MPI.Waitall!(reqs)

        # Unload the buffer.
        for i in 1:p.npx
            js = (i-1)*g.jblock + 1; je = i*g.jblock
            recvbuf_view = @view recvbuf[:, :, :, i]
            @tturbo data_out[:, js:je, :] .= recvbuf_view
        end
    end
end


function transpose_yx!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelDistributed)
    if data_out === data
        return
    elseif size(data_out) == size(data)
        @tturbo data_out[:, :, :] .= data
    else
        sendbuf = reshape(sendbuf, (g.iblock, g.jmax, g.kblock, p.npy))
        recvbuf = reshape(recvbuf, (g.iblock, g.jmax, g.kblock, p.npy))

        # Load the buffer.
        for i in 1:p.npy
            js = (i-1)*g.jmax + 1; je = i*g.jmax
            data_view = @view data[:, js:je, :]
            @tturbo sendbuf[:, :, :, i] .= data_view
        end

        # Communicate data.
        # message_size = g.iblock*g.jmax*g.kblock
        # MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commy)

        # Communicate data using point-to-point.
        reqs = Vector{MPI.Request}(undef, 2*p.npy)

        for i in 1:p.npy
            recvbuf_view = @view recvbuf[:, :, :, i];
            reqs[i] = MPI.Irecv!(recvbuf_view, i-1, 1, p.commy)
        end

        for i in 1:p.npy
            sendbuf_view = @view sendbuf[:, :, :, i]
            reqs[i+p.npy] = MPI.Isend(sendbuf_view, i-1, 1, p.commy)
        end

        MPI.Waitall!(reqs)

        # Unload the buffer.
        for i in 1:p.npy
            is = (i-1)*g.iblock + 1; ie = i*g.iblock
            recvbuf_view = @view recvbuf[:, :, :, i]
            @tturbo data_out[is:ie, :, :] .= recvbuf_view
        end
    end
end


function transpose_xz!(data_out, data, sendbuf, recvbuf, g::Grid, p::ParallelDistributed)
    if data_out === data
        return
    elseif size(data_out) == size(data)
        @tturbo data_out[:, :, :] .= data
    else
        sendbuf = reshape(sendbuf, (g.imax, g.jmax, g.kblock, p.npx))
        recvbuf = reshape(recvbuf, (g.imax, g.jmax, g.kblock, p.npx))

        # Load the buffer.
        for i in 1:p.npx
            is = (i-1)*g.imax + 1; ie = i*g.imax
            data_view = @view data[is:ie, :, :]
            @tturbo sendbuf[:, :, :, i] .= data_view
        end

        # Communicate data.
        # message_size = g.imax*g.jmax*g.kblock
        # MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commx)

        # Communicate data using point-to-point.
        reqs = Vector{MPI.Request}(undef, 2*p.npx)

        for i in 1:p.npx
            recvbuf_view = @view recvbuf[:, :, :, i];
            reqs[i] = MPI.Irecv!(recvbuf_view, i-1, 1, p.commx)
        end

        for i in 1:p.npx
            sendbuf_view = @view sendbuf[:, :, :, i]
            reqs[i+p.npx] = MPI.Isend(sendbuf_view, i-1, 1, p.commx)
        end

        MPI.Waitall!(reqs)

        # Unload the buffer.
        for i in 1:p.npx
            ks = (i-1)*g.kblock + 1; ke = i*g.kblock
            recvbuf_view = @view recvbuf[:, :, :, i]
            @tturbo data_out[:, :, ks:ke] .= recvbuf_view
        end
    end
end

