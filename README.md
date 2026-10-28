# letter-canonical

Code and figure scripts for a letter on guiding-center orbit integration in
global canonical coordinates, comparing symplectic and non-symplectic
integrators in tokamak geometry.

The manuscript is kept on Overleaf; `writeup_canonical.lyx` and
`canonical_comment.lyx` hold accompanying notes, and `canonical.bib` the
bibliography.

## Contents

- `SRC/`: Fortran sources for magnetic field evaluation (test field,
  cylindrical and tokamak equilibria, Biot-Savart coils), canonical field
  representations and orbit integrators (RK45, DOP853, explicit-implicit Euler,
  symplectic midpoint)
- `RUN/`: input files and Python scripts for the letter's figures
- `PLOT/`: further plotting scripts
- `DATA/`: equilibrium (g-file) and coil/wall data used by the runs
- `test/`, `SRC/test_*.f90`: tests
- `gen/`: generator for `SRC/generated/segment_log.inc` (see
  [gen/README.md](gen/README.md))

## Building

Requirements: CMake 3.18 or newer, Ninja, a Fortran compiler (GNU, Intel or
NVHPC) with OpenMP, BLAS and LAPACK, and Python 3 for the tests. Selected
[libneo](https://github.com/itpplasma/libneo) sources are fetched at a pinned
commit during configuration.

    make            # release build in BUILD/
    make debug
    cd BUILD && ctest --output-on-failure

## Reproducing the figures

    cd RUN
    make            # runs BUILD/letter_canonical.x on the fig*_*.in cases
    python fig1_momentum_energy.py
    python fig2_islands.py
    python fig3_bananatips.py

`RUN/Makefile` lists the cases it runs; further `fig*_*.in` inputs can be run
with `./run.sh <input>.in`.

- Fig. 1: banana orbit crossing the separatrix in axisymmetric geometry.
  Symplectic schemes conserve toroidal momentum exactly and energy within
  bounds. All integrators used about 2000 field evaluations per bounce period;
  explicit-implicit Euler could increase the step length further, reducing
  field evaluations by a factor 3.9, before reaching its stability limit.
- Fig. 2: islands of passing orbits.
- Fig. 3: banana tips.

## Paper

    make paper

clones the Overleaf project into `paper/` (or pulls it if present). Log in with
`git` as user name and an Overleaf authentication token as password; access to
the project is required.

## Output format

Callback output files contain time and the five phase-space coordinates in the
first six columns and the orbit index in the seventh. Records from parallel
runs may be interleaved, so group by orbit index before plotting a trajectory.

Cylindrical field evaluations serialise access to libneo's shared state, so the
OpenMP loop over orbits is thread-safe but not necessarily faster. The
`parallel_orbits` test compares the callback records of eight synthetic orbits
run with one and four threads.

## License

MIT, see [LICENSE](LICENSE).
