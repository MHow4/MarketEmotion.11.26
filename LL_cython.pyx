# cython: language_level=3
import numpy as np
cimport numpy as cnp

cnp.import_array()


def one_energy_checkerboard(cnp.ndarray centre, tuple nbrs):
    # centre: array of angles at the sites; nbrs: tuple of 4 neighbour arrays
    cdef cnp.ndarray en = np.zeros_like(centre)
    cdef cnp.ndarray nb
    for nb in nbrs:
        en = en + 0.5 * (1.0 - 3.0 * np.cos(centre - nb) ** 2)
    return en


def MC_step(cnp.ndarray arr, double Ts, int nmax):
    # constants
    cdef double scale = 0.1 + Ts
    cdef long accept = 0

    # create checkerboard masks
    cdef cnp.ndarray i, j
    i, j = np.indices((nmax, nmax))
    cdef list masks = [(i + j) % 2 == 0,
                       (i + j) % 2 == 1]

    cdef cnp.ndarray mask, aran, en0, en1, dE, rand, ok
    cdef tuple nbrs

    for mask in masks:
        aran = np.random.normal(scale=scale, size=arr.shape)
        nbrs = (np.roll(arr, 1, axis=0), np.roll(arr, -1, axis=0),
                np.roll(arr, 1, axis=1), np.roll(arr, -1, axis=1))

        en0 = one_energy_checkerboard(arr, nbrs)
        en1 = one_energy_checkerboard(arr + aran, nbrs)
        dE = en1 - en0

        # Metropolis test: accept if dE <= 0, else with prob exp(-dE/Ts)
        rand = np.random.uniform(0.0, 1.0, size=arr.shape)
        ok = (dE <= 0) | (np.exp(-np.maximum(dE, 0) / Ts) >= rand)
        # only update the current sublattice
        ok &= mask

        # update accepted sites
        arr += np.where(ok, aran, 0.0)
        accept += ok.sum()

    return accept / (nmax * nmax)
