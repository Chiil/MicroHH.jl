abstract type Parallel end


struct ParallelDistributed <: Parallel
    npx::Int64
    npy::Int64

    id::Int64
    idx::Int64
    idy::Int64

    commxy
    commx
    commy
end


struct ParallelSerial <: Parallel
    npx::Int64
    npy::Int64

    id::Int64
    idx::Int64
    idy::Int64

    commxy
    commx
    commy
end


function Parallel(npx, npy)
    if npx > 1 || npy > 1
        MPI.Init()
        
        dims = [npy, npx]; periodic = [1, 1]; reorder = true

        commxy = MPI.Cart_create(MPI.COMM_WORLD, dims, periodic, reorder)
        commx = MPI.Cart_sub(commxy, [false, true])
        commy = MPI.Cart_sub(commxy, [true, false])

        id = MPI.Comm_rank(commxy)
        idx = MPI.Comm_rank(commx)
        idy = MPI.Comm_rank(commy)

        return ParallelSerial(npx, npy, id, idx, idy, commxy, commx, commy)
    else
        id = 0; idx = 0; idy = 0
        commxy = Nothing; commx = Nothing; commy = Nothing
    end

    return ParallelDistributed(npx, npy, id, idx, idy, commxy, commx, commy)
end


function close_parallel(p::ParallelSerial)
end


function close_parallel(p::ParallelDistributed)
    MPI.Finalize()
end
