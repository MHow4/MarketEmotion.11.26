# Lebwohl-Lasher Monte Carlo Simulation (2D)

The simulation runs on a 2D square lattice. Each cell holds a single angle describing the orientation of a rod-like (nematic) molecule. Neighbouring cells interact through a second-order Legendre-polynomial potential, and the system is evolved with the Metropolis algorithm.

## Branches

| Branch    | Description                                                                    |
|-----------|--------------------------------------------------------------------------------|
| `cython`  | Uses the vectorised code and cython to translate the code into C               |
| `master`  | MAIN Branch - Contains the original Lebwohl Lasher simulation                  |
| `mpi4py`  | Uses the vectorised code, along with MPI to parallelise the algorithm          |
| `numba`   | Uses just in time compilation with Numba on the original code                  |
| `vector`  | Uses NumPy to calculate the MC step across the array in a checkerboard pattern |

To switch between branches, use:

```bash
git branch <BRANCH>
```

## Usage

### Master, Numba, Vector

```bash
python LebwohlLasher.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG>
```

| Argument      | Type  | Description                                                                                          |
|---------------|-------|------------------------------------------------------------------------------------------------------|
| `ITERATIONS`  | int   | Number of Monte Carlo steps (MCS). 1 MCS = `SIZE*SIZE` attempted moves, i.e. one attempt per cell on average. |
| `SIZE`        | int   | Side length of the square lattice (the lattice is `SIZE x SIZE`).                                    |
| `TEMPERATURE` | float | Reduced temperature `T* = kT/ε`, in the range 0.0 - 2.0.                                             |
| `PLOTFLAG`    | int   | `0` = no plot, `1` = energy plot, `2` = angle plot, any other value = plain black plot.              |

### Example

```bash
python LebwohlLasher.py 100 50 0.5 0
```

Runs 100 MC steps on a 50x50 lattice at `T* = 0.5` with no plotting.

Example console output:

```
LebwohlLasher.py: Size: 50, Steps: 100, T*: 0.500: Order: 0.xxx, Time: x.xxxxxx s
```

### MPI4Py

```bash
mpiexec -n <THREADS> python LebwohlLasher.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG>
```

To run the MPI code, add the ```mpiexec -n <THREADS>``` prefix to the original running script. ```<THREADS>``` determines the number of threads used

### Cython

In order to run the Cython file, first the Cython needs to be compiled using:

```bash
python setup.py build_ext --inplace
```

Once this has been done, the code can be run using:

```bash
$env:OMP_NUM_THREADS = <THREADS>
python LebwohlLasher.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG>
```

## Output

### Console

A one-line summary with lattice size, number of steps, temperature, order parameter and run time.

### Data file

A text file named automatically from the date and time at the end of the run:

```
LL-Output-<Day>-<DD>-<Mon>-<YYYY>-at-<HH>-<MM>-<SS><AM/PM>.txt
```

It contains a header with the run parameters (lattice size, number of steps, reduced temperature, run time) followed by one row per MC step (including step 0, the initial state) with the columns:

| Column    | Description                                      |
|-----------|--------------------------------------------------|
| MC step   | Step number                                      |
| Ratio     | Acceptance ratio (fraction of accepted moves)    |
| Energy    | Total reduced energy of the lattice              |
| Order     | Order parameter `S`                              |

The acceptance ratio for step 0 is set to a nominal value of 0.5. As a rule of thumb, aim to keep the acceptance ratio around 0.5 for an efficient simulation.


## Reference

P.A. Lebwohl and G. Lasher, "Nematic-Liquid-Crystal Order - A Monte Carlo Calculation," *Phys. Rev. A* **6**, 426 (1972).

Original code: SH, 16-Oct-23.
