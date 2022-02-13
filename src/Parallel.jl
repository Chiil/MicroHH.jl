abstract type Parallel end


struct ParallelDistributed <: Parallel
    npx::Int64
    npy::Int64

    id::Int64

    id_x::Int64
    id_y::Int64

    id_west::Int64
    id_east::Int64
    id_south::Int64
    id_north::Int64

    commxy
    commx
    commy
end


struct ParallelSerial <: Parallel
    npx::Int64
    npy::Int64

    id::Int64

    id_x::Int64
    id_y::Int64

    id_west::Int64
    id_east::Int64
    id_south::Int64
    id_north::Int64

    commxy
    commx
    commy
end


function Parallel(npx, npy)
    if do_mpi[]
        MPI.Init()

        dims = [npy, npx]; periodic = [1, 1]; reorder = true

        commxy = MPI.Cart_create(MPI.COMM_WORLD, dims, periodic, reorder)
        commx = MPI.Cart_sub(commxy, [false, true])
        commy = MPI.Cart_sub(commxy, [true, false])

        id = MPI.Comm_rank(commxy)
        id_x = MPI.Comm_rank(commx)
        id_y = MPI.Comm_rank(commy)

        id_west, id_east = MPI.Cart_shift(commxy, 1, 1)
        id_south, id_north = MPI.Cart_shift(commxy, 0, 1)

        if id == 0
            @info "Initialized MicroHH with do_mpi = $(do_mpi[]), on npx = $(npx[]), npy = $(npy[]) tasks, nthreads = $(Threads.nthreads())."
        end

        return ParallelDistributed(
            npx, npy,
            id, id_x, id_y,
            id_west, id_east, id_south, id_north,
            commxy, commx, commy)
    else
        id = 0; id_x = 0; id_y = 0
        id_west = 0; id_east = 0; id_south = 0; id_north = 0
        commxy = Nothing; commx = Nothing; commy = Nothing

        @info "Initialized MicroHH with do_mpi = $(do_mpi[]), on npx = $(npx[]), npy = $(npy[]) tasks, nthreads = $(Threads.nthreads())."

        return ParallelSerial(
            npx, npy,
            id, id_x, id_y,
            id_west, id_east, id_south, id_north,
            commxy, commx, commy)
    end
end


function close_parallel(p::ParallelSerial)
end


function close_parallel(p::ParallelDistributed)
    MPI.Finalize()
end
