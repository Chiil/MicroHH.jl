## User input.
npx = 2; npy = 4
itot = 256; jtot = 128; ktot = 64

imax = itot ÷ npx; jmax = jtot ÷ npy


## Init MPI and create grid.
using MPI

MPI.Init()

dims = [npx, npy]; periodic = [1, 1]; reorder = true
commxy = MPI.Cart_create(MPI.COMM_WORLD, dims, periodic, reorder)
mpiid = MPI.Comm_rank(commxy)

nwest, neast = MPI.Cart_shift(commxy, 1, 1)
nsouth, nnorth = MPI.Cart_shift(commxy, 0, 1)

data = ones(Int, imax+2, jmax+2, ktot+2) * mpiid

function boundary_cyclic(data)
    # Transfer the east-west data.
    send_east = data[end-1, :, :]
    send_west = data[2, :, :]
    recv_west = zeros(Int, size(send_east))
    recv_east = zeros(Int, size(send_west))

    reqs = Vector{MPI.Request}(undef, 4)
    reqs[1] = MPI.Irecv!(recv_west, nwest, 1, commxy);
    reqs[2] = MPI.Irecv!(recv_east, neast, 2, commxy);
    reqs[3] = MPI.Isend( send_east, neast, 1, commxy);
    reqs[4] = MPI.Isend( send_west, nwest, 2, commxy);
    MPI.Waitall!(reqs)

    data[1, :, :] .= recv_west[:, :]
    data[end, :, :] .= recv_east[:, :]

    # Transfer the north-south data.
    send_north = data[:, end-1, :]
    send_south = data[:, 2, :]
    recv_south = zeros(Int, size(send_north))
    recv_north = zeros(Int, size(send_south))

    reqs[1] = MPI.Irecv!(recv_north, nnorth, 1, commxy);
    reqs[2] = MPI.Irecv!(recv_south, nsouth, 2, commxy);
    reqs[3] = MPI.Isend( send_south, nsouth, 1, commxy);
    reqs[4] = MPI.Isend( send_north, nnorth, 2, commxy);
    MPI.Waitall!(reqs)

    data[:, 1, :] .= recv_south[:, :]
    data[:, end, :] .= recv_north[:, :]
end


## Run benchmark.
boundary_cyclic(data)

for n in 1:10
    dt = @elapsed boundary_cyclic(data)
    if mpiid == 0 println("Elapsed: $dt (s)") end
end


# Print output.
if mpiid == 4
    print("Neighbors (W, E): $mpiid = $nwest, $neast\n")
    print("Neighbors (S, N): $mpiid = $nsouth, $nnorth\n")
    print("$(data[:, 1, 1])\n")
end
MPI.Barrier(commxy)
