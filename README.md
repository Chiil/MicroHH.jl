# MicroHH

This is a Julia port of the dynamical core of [MicroHH](https://microhh.github.io), a 3D CFD code to simulate the atmospheric boundary layer. It is made for educational and exploratory research purposes, including machine learning. In comparison to its C++/CUDA reference it lacks all physical parameterizations, but it comes with the ability to run multiple simulations at the same time and let them interact.

A few example cases are in the `cases` folder.

MicroHH should be launched without flags for serial mode:

* `drycbl` with 1 thread(s): `julia --project -O3 drycbl_run.nl`
* `drycbl` with 2 thread(s): `julia --project -O3 -t2 drycbl_run.nl`

MicroHH should be launched with `--use-mpi` flags for parallel mode with `--npx` and `--npy` specifying the MPI tiling:

* `drycbl` with 8 tasks and 1 thread(s): `mpiexecjl -n 8 julia --project -O3 drycbl_run.nl --use-mpi --npx 2 --npy 4`
* `drycbl` with 4 tasks and 2 thread(s): `mpiexecjl -n 4 julia --project -O3 -t2 drycbl_run.nl --use-mpi --npx 2 --npy 2`
* `drycbl` with 2 tasks and 4 thread(s): `mpiexecjl -n 2 julia --project -O3 -t4 drycbl_run.nl --use-mpi --npx 1 --npy 2`

