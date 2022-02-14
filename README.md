# MicroHH

This is a Julia port of the dynamical core of [MicroHH](https://microhh.github.io), a 3D CFD code to simulate the atmospheric boundary layer. It is made for educational and exploratory research purposes, including machine learning. In comparison to its C++/CUDA reference it lacks all physical parameterizations, but it comes with the ability to run multiple simulations at the same time and let them interact.

A few example cases are in the `cases` folder.

MicroHH should be launched without flags for serial mode.

Example for `drycbl` case with 1 threads: `julia --project -O3 drycbl_run.nl`
Example for `drycbl` case with 2 threads: `julia --project -O3 -t2 drycbl_run.nl`

MicroHH should be launched with `--use-mpi` flags for parallel mode with `--npx` and `--npy` specifying the MPI tiling.

Example for `drycbl` case with 8 tasks and 1 threads: `mpiexecjl -n 4 julia --project -O3 drycbl_run.nl --use-mpi --npx 2 --npy 4`
Example for `drycbl` case with 4 tasks and 2 threads: `mpiexecjl -n 4 julia --project -O3 -t2 drycbl_run.nl --use-mpi --npx 2 --npy 2`

