# MicroHH

This is a Julia port of the dynamical core of [MicroHH](https://microhh.github.io), a 3D CFD code to simulate the atmospheric boundary layer. It is made for educational and exploratory research purposes, including machine learning. In comparison to its C++/CUDA reference it lacks all physical parameterizations, but it comes with the ability to run multiple simulations at the same time and let them interact.

A few example cases are in the `cases` folder.

MicroHH is by default precompiled for serial computing with Float64 precision. It can be run on a chosen number of threads:

* `drycbl` with 1 thread(s): `julia --project -O3 drycbl_run.nl`
* `drycbl` with 2 thread(s): `julia --project -O3 -t2 drycbl_run.nl`

If desired, MPI mode can be enabled via de `set_use_mpi!` function that can be run from Julia, and the settings are stored in `LocalPreferences.toml`.

```
using MicroHH
MicroHH.set_use_mpi!(true)
```

With MPI, the simulation should be started for instance as (if `npx = 2` and `npy=2`)
`drycbl` with 4 tasks and 2 thread(s): `mpiexec -n 4 julia --project -O3 -t2 drycbl_run.nl`

MicroHH can be run at single precision by setting the `use_sp` flag to `true`.

```
using MicroHH
MicroHH.set_use_sp!(true)
```
