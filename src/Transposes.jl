function transpose_zx(data_out, data, g::Grid, p::ParallelDistributed)
    sendbuf = reshape(similar(data), (g.imax, g.jmax, g.kblock, p.npx))
    recvbuf = similar(sendbuf)

    # Load the buffer.
    for i in 1:p.npx
        ks = (i-1)*g.kblock + 1; ke = i*g.kblock
        @tturbo sendbuf[:, :, :, i] .= data[:, :, ks:ke]
    end

    # Communicate data.
    message_size = g.imax*g.jmax*g.kblock
    MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commx)

    # Unload the buffer.
    for i in 1:p.npx
        is = (i-1)*g.imax + 1; ie = i*g.imax
        @tturbo data_out[is:ie, :, :] .= recvbuf[:, :, :, i]
    end
end


function transpose_xy(data_out, data, g::Grid, p::ParallelDistributed)
    sendbuf = reshape(similar(data), (g.iblock, g.jmax, g.kblock, p.npy))
    recvbuf = similar(sendbuf)

    # Load the buffer.
    for i in 1:p.npy
        is = (i-1)*g.iblock + 1; ie = i*g.iblock
        @tturbo sendbuf[:, :, :, i] .= data[is:ie, :, :]
    end

    # Communicate data.
    message_size = g.iblock*g.jmax*g.kblock
    MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commy)

    # Unload the buffer.
    for i in 1:p.npy
        js = (i-1)*g.jmax + 1; je = i*g.jmax
        @tturbo data_out[:, js:je, :] .= recvbuf[:, :, :, i]
    end
end


function transpose_yzt(data_out, data, g::Grid, p::ParallelDistributed)
    sendbuf = reshape(similar(data), (g.iblock, g.jblock, g.kblock, p.npx))
    recvbuf = similar(sendbuf)

    # Load the buffer.
    for i in 1:p.npx
        js = (i-1)*g.jblock + 1; je = i*g.jblock
        @tturbo sendbuf[:, :, :, i] .= data[:, js:je, :]
    end

    # Communicate data.
    message_size = g.iblock*g.jblock*g.kblock
    MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commx)

    # Unload the buffer.
    for i in 1:p.npx
        ks = (i-1)*g.kblock + 1; ke = i*g.kblock
        @tturbo data_out[:, :, ks:ke] .= recvbuf[:, :, :, i]
    end
end


function transpose_zty(data_out, data, g::Grid, p::ParallelDistributed)
    sendbuf = reshape(similar(data), (g.iblock, g.jblock, g.kblock, p.npx))
    recvbuf = similar(sendbuf)

    # Load the buffer.
    for i in 1:p.npx
        ks = (i-1)*g.kblock + 1; ke = i*g.kblock
        @tturbo sendbuf[:, :, :, i] .= data[:, :, ks:ke]
    end

    # Communicate data.
    message_size = g.iblock*g.jblock*g.kblock
    MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commx)

    # Unload the buffer.
    for i in 1:p.npx
        js = (i-1)*g.jblock + 1; je = i*g.jblock
        @tturbo data_out[:, js:je, :] = recvbuf[:, :, :, i]
    end
end


function transpose_yx(data_out, data, g::Grid, p::ParallelDistributed)
    sendbuf = reshape(similar(data), (g.iblock, g.jmax, g.kblock, p.npy))
    recvbuf = similar(sendbuf)

    # Load the buffer.
    for i in 1:p.npy
        js = (i-1)*g.jmax + 1; je = i*g.jmax
        @tturbo sendbuf[:, :, :, i] .= data[:, js:je, :]
    end

    # Communicate data.
    message_size = g.iblock*g.jmax*g.kblock
    MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commy)

    # Unload the buffer.
    for i in 1:p.npy
        is = (i-1)*g.iblock + 1; ie = i*g.iblock
        @tturbo data_out[is:ie, :, :] .= recvbuf[:, :, :, i]
    end
end


function transpose_xz(data_out, data, g::Grid, p::ParallelDistributed)
    sendbuf = reshape(similar(data), (g.imax, g.jmax, g.kblock, p.npx))
    recvbuf = similar(sendbuf)

    # Load the buffer.
    for i in 1:p.npx
        is = (i-1)*g.imax + 1; ie = i*g.imax
        @tturbo sendbuf[:, :, :, i] .= data[is:ie, :, :]
    end

    # Communicate data.
    message_size = g.imax*g.jmax*g.kblock
    MPI.Alltoall!(MPI.UBuffer(sendbuf, message_size), MPI.UBuffer(recvbuf, message_size), p.commx)

    # Unload the buffer.
    for i in 1:p.npx
        ks = (i-1)*g.kblock + 1; ke = i*g.kblock
        @tturbo data_out[:, :, ks:ke] .= recvbuf[:, :, :, i]
    end
end
