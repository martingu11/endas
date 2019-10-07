"""
Cython implementation of endas/localization/taper.py
"""
import cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def taper_gc(double[:] x not None, double[:] d not None, double[:] out not None, double L):
    # At this point we assume all the necessary asserts were already done!
    cdef int n = x.shape[0]
    cdef int i
    cdef double r
    for i in range(n):
        r = d[i] / L
        if r < 1: out[i] = x[i] * 1.0 - (5/3.0)*r**2 + (5/8.0)*r**3 + (1/2.0)*r**4 - (1/4.0)*r**5
        elif r < 2: out[i] = x[i] * 4.0 - 5.0*r + (5/3.0)*r**2 + (5/8.0)*r**3 - (1/2.0)*r**4 + (1/12.0)*r**5 - (2.0/3.0)*r
        else: out[i] = 0



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def taper_linear(double[:] x not None, double[:] d not None, double[:] out not None, double L):
    # At this point we assume all the necessary asserts were already done!

    cdef int n = x.shape[0]
    cdef int i
    cdef double r
    for i in range(n):
        r = d[i] / L
        if r < 1.0: out[i] = x[i] * 1.0 - r
        else: out[i] = 0



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def taper_spherical(double[:] x not None, double[:] d not None, double[:] out not None, double L):
    # At this point we assume all the necessary asserts were already done!
    cdef int n = x.shape[0]
    cdef int i
    cdef double r
    for i in range(n):
        r = d[i] / L
        if r < 1.0: out[i] = x[i] * 1.0 - (1.5*r - 0.5*r**3)
        else: out[i] = 0



