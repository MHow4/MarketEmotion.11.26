# cython: language_level=3, boundscheck=False, wraparound=False
# cython: cdivision=True, initializedcheck=False
import numpy as np
from libc.math cimport cos, sin, exp
from cython.parallel cimport prange

cdef inline double pair_energy(double a, double b) noexcept nogil:
    cdef double c = cos(a - b)
    return 0.5 * (1.0 - 3.0 * c * c)

cdef double one_energy(double[:, ::1] arr, int ix, int iy, int nmax) noexcept nogil:
    cdef int ixp = ix + 1 if ix < nmax - 1 else 0
    cdef int ixm = ix - 1 if ix > 0 else nmax - 1
    cdef int iyp = iy + 1 if iy < nmax - 1 else 0
    cdef int iym = iy - 1 if iy > 0 else nmax - 1

    cdef double c = arr[ix, iy]
    return (pair_energy(c, arr[ixp, iy]) + pair_energy(c, arr[ixm, iy]) +
            pair_energy(c, arr[ix, iyp]) + pair_energy(c, arr[ix, iym]))


def MC_step(double[:, ::1] arr, double Ts, int nmax):
    # constants
    cdef double scale = 0.1 + Ts
    cdef long accept = 0

    cdef int parity, i, j, j0
    cdef double old, e0, e1, dE
    cdef double[:, ::1] aran, rnd

    for parity in range(2):
        # random numbers are generated up front
        aran = np.random.normal(scale=scale, size=(nmax, nmax))
        rnd = np.random.uniform(0.0, 1.0, size=(nmax, nmax))

        # sites with (i+j)%2 == parity have no neighbours so can be updated in parallel.
        for i in prange(nmax, nogil=True, schedule='static'):

            j0 = (i + parity) % 2
            for j in range(j0, nmax, 2):

                # calculate energies
                old = arr[i, j]
                e0 = one_energy(arr, i, j, nmax)
                arr[i, j] = old + aran[i, j]
                e1 = one_energy(arr, i, j, nmax)
                dE = e1 - e0

                # energy update check
                if dE <= 0.0 or exp(-dE / Ts) >= rnd[i, j]:
                    accept += 1
                else:
                    arr[i, j] = old

    return accept / (<double>nmax * nmax)
