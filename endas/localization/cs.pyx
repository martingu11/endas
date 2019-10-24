"""
Coordinate Systems.
"""

__all__ = ['EuclideanCS', 'LatLonCS']

import numpy as np
import cython
from libc.math cimport sqrt, sin, cos, asin

from . import CoordinateSystem


ctypedef fused coord_t:
    int
    float
    double


# Euclidean distance specialized for the most common cases - 1d, 2d and 3d.
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_euclid_distance_1d(coord_t[:] A not None, coord_t[:] B not None, coord_t[:] out not None):
    # At this point we assume all the necessary asserts were already done! A and B are flat arrays
    cdef int n = A.shape[0]
    cdef int i
    cdef coord_t d
    for i in range(n):
        d = A[i] - B[i]
        if d < 0: d = -d   # Does anyone know better way to do abs() for both ints and floats?
        out[i] = d


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_euclid_distance_2d(double[:,:] A not None, double[:,:] B not None, double[:] out not None):
    # At this point we assume all the necessary asserts were already done! A and B have shape (n, 2)
    cdef int n = A.shape[0]
    cdef int i
    for i in range(n):
        out[i] = sqrt( (A[i,0]-B[i,0])**2 + (A[i,1]-B[i,1])**2 )


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_euclid_distance_3d(coord_t[:,:] A not None, coord_t[:,:] B not None, double[:] out not None):
    # At this point we assume all the necessary asserts were already done! A and B have shape (n, 3)
    cdef int n = A.shape[0]
    cdef int i
    for i in range(n):
        out[i] = sqrt( (A[i,0]-B[i,0])**2 + (A[i,1]-B[i,1])**2 + (A[i,2]-B[i,2])**2 )


# Euclidean distance for generic N-d case with extra inner loop
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_euclid_distance_Nd(coord_t[:,:] A not None, coord_t[:,:] B not None, double[:] out not None):
    # At this point we assume all the necessary asserts were already done! A and B have shape (n, N)
    cdef int n = A.shape[0]
    cdef int N = A.shape[1]
    cdef int i
    cdef double s
    for i in range(n):
        s = 0.0
        for j in range(N): s+= (A[i,j] - B[i,j])**2
        out[i] = sqrt(s)


class EuclideanCS(CoordinateSystem):
    """
    Cartesian coordinate system in N-dimensional Euclidean space.

    Args:
        ndim : Dimensionality of the coordinate system, must be >= 1

    """

    def __init__(self, ndim, dtype=np.double):
        assert ndim >= 1
        self._ndim = ndim
        self._dtype = dtype

    @property
    def ndim(self): return self._ndim

    @property
    def is_cartesian(self): return True

    def distance(self, A, B, out=None):
        assert A.shape[1] == self.ndim
        assert A.shape == B.shape or A.shape[0] <= 1

        n = A.shape[0]
        if n == 0: return np.empty(0)

        if out is None: out = np.empty(n, dtype=self._dtype)
        if self.ndim == 1:   cs_euclid_distance_1d(A, B, out)
        elif self.ndim == 2: cs_euclid_distance_2d(A, B, out)
        elif self.ndim == 3: cs_euclid_distance_3d(A, B, out)
        else:                cs_euclid_distance_Nd(A, B, out)





# Lat-lon distance on a perfect sphere using the Haversine formula
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_latlon_distance(double[:,:] A not None, double[:,:] B not None, double[:] out not None, double R):
    # At this point we assume all the necessary asserts were already done! A and B have shape (n, 2)
    cdef int n = A.shape[0]
    cdef int i
    cdef double Alat, Blat, Alon, Blon
    cdef double dfif_lat, diff_lon
    cdef double a

    for i in range(n):
        Alat = A[i,0] * 0.0174532925    # Degree to Radian
        Alon = A[i,1] * 0.0174532925
        Blat = B[i,0] * 0.0174532925
        Blon = B[i,1] * 0.0174532925

        a = sin((Alat - Blat) / 2.0)**2 + cos(Alat) * cos(Blat) * sin((Alon - Blon) / 2.0)**2
        if a > 1.0: a = 1.0
        a = 2.0 * asin(sqrt(a))

        out[i] = a * R



class LatLonCS(CoordinateSystem):
    """
    Polar coordinate system on a perfect sphere.

    ``LatLonCS`` implementes coordinate system on a perfect sphere with coordinates of any point expressed as latitude
    and longitude. The currently implemented coordinate system is strictly two-dimensional, i.e. there is no vertical
    component. The coordinates are assumed to be in degrees.

    Args:
        R : The radius of the great circle on which the distance is calculated. The default value of 6371 km corresponds
            to the mean Earth radius (R1) as defined by the International Union of Geodesy and Geophysics.

    This class is a simple implementation of a polar coordinate system assuming a perfect sphere of radius *R* and uses
    the `Haversine formula <https://en.wikipedia.org/wiki/Haversine_formula>`_ to compute the great-circle distance.
    If the use of a better approximation is required (i.e. a spheroid), please consider coding your own.
    """
    def __init__(self, R=6.371e6):
        self.R = R

    @property
    def ndim(self): return 2

    @property
    def is_cartesian(self): return False

    def distance(self, A, B, out=None):
        assert A.shape == B.shape
        assert A.shape[1] == self.ndim
        n = A.shape[0]
        if n == 0: return np.empty(0)

        if out is None: out = np.empty(n, dtype=np.double)
        cs_latlon_distance(A, B, out, self.R)
        return out



