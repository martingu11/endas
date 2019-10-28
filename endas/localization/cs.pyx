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



class EuclideanCS(CoordinateSystem):
    """
    Cartesian coordinate system in N-dimensional Euclidean space.

    Args:
        ndim : Dimensionality of the coordinate system, must be >= 1

    """

    def __init__(self, ndim):
        assert ndim >= 1
        self._ndim = ndim


        #if self.ndim == 1: self._fn = self._distance_1d
        #elif self.ndim == 2: self._fn = self._distance_2d
        #elif self.ndim == 3: self._fn = self._distance_3d
        #elif self.ndim == 4: self._fn = self._distance_Nd

    @property
    def ndim(self): return self._ndim

    @property
    def is_cartesian(self): return True

    def distance(self, A, B, out=None):

        have_single_A = A.shape[0] == 1
        if have_single_A: assert A.shape[1] == self.ndim
        else: assert A.shape == B.shape
        assert B.shape[1] == self.ndim

        n = B.shape[0]
        if n == 0: return np.empty(0)

        if out is None:
            out = np.empty(n, dtype=np.double)
        else:
            assert out.ndim == 1
            assert out.size == n

        if self.ndim == 1: self._distance_1d(A.ravel(), B.ravel(), out)
        elif self.ndim == 2: self._distance_2d(A, B, out)
        elif self.ndim == 3: self._distance_3d(A, B, out)
        else: self._distance_Nd(A, B, out)

        return out


    # One-dimensional case
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def _distance_1d(self, coord_t[::1] A, coord_t[::1] B, double[::1] out):
        cdef int n = B.shape[0]
        cdef int i
        if A.shape[0] == 1:
            for i in range(n):
                out[i] = abs(A[0] - B[i])
        else:
            for i in range(n):
                out[i] = abs(A[i] - B[i])


    # Two-dimensional case
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def _distance_2d(self, coord_t[:,::1] A not None, coord_t[:,::1] B not None, double[::1] out not None):
        cdef int n = B.shape[0]
        cdef int i
        if A.shape[0] == 1:
            for i in range(n):
                out[i] = sqrt( (A[0,0] - B[i,0])**2 + (A[0,1] - B[i,1])**2 )
        else:
            for i in range(n):
                out[i] = sqrt( (A[i,0] - B[i,0])**2 + (A[i,1] - B[i,1])**2 )

    # Three-dimensional case
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def _distance_3d(self, coord_t[:,::1] A not None, coord_t[:,::1] B not None, double[::1] out not None):
        cdef int n = B.shape[0]
        cdef int i
        if A.shape[0] == 1:
            for i in range(n):
                out[i] = sqrt( (A[0,0] - B[i,0])**2 + (A[0,1] - B[i,1])**2 + (A[0,2] - B[i,2])**2 )
        else:
            for i in range(n):
                out[i] = sqrt( (A[i,0] - B[i,0])**2 + (A[i,1] - B[i,1])**2 + (A[i,2] - B[i,2])**2)


    # Generic N-dimensional case
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def _distance_Nd(self, coord_t[:,::1] A not None, coord_t[:,::1] B not None, double[::1] out not None):
        cdef int n = B.shape[0]
        cdef int N = B.shape[1]
        cdef int i
        cdef double sum_sq
        if A.shape[0] == 1:
            for i in range(n):
                sum_sq = 0.0
                for j in range(N): sum_sq += (A[0,j] - B[i,j])**2
                out[i] = sqrt(sum_sq)
        else:
            for i in range(n):
                sum_sq = 0.0
                for j in range(N): sum_sq += (A[i,j] - B[i,j])**2
                out[i] = sqrt(sum_sq)




# Haversine formula for great circle distance on sphere
@cython.cdivision(True)
cdef inline double haversine(double Alat, double Alon, double Blat, double Blon, double R) nogil:
    Alat*= 0.0174532925 # Conversion to radians
    Blat*= 0.0174532925
    Alon*= 0.0174532925
    Blon*= 0.0174532925

    cdef double a = sin((Alat - Blat) / 2.0)**2 + cos(Alat) * cos(Blat) * sin((Alon - Blon) / 2.0)**2
    if a > 1.0: a = 1.0
    a = 2.0 * asin(sqrt(a))
    return a * R


class LatLonCS(CoordinateSystem):
    """
    Polar coordinate system on a perfect sphere.

    ``LatLonCS`` implements coordinate system on a perfect sphere with coordinates of any point expressed as latitude
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
        have_single_A = A.shape[0] == 1
        if have_single_A: assert A.shape[1] == self.ndim
        else: assert A.shape == B.shape
        assert B.shape[1] == self.ndim

        n = B.shape[0]
        if n == 0: return np.empty(0)

        if out is None:
            out = np.empty(n, dtype=np.double)
        else:
            assert out.ndim == 1
            assert out.size == n

        cdef int i
        cdef double Alat, Blat, Alon, Blon
        cdef double dfif_lat, diff_lon
        cdef double a

        cdef double[:,::1] A_view = A
        cdef double[:,::1] B_view = B
        cdef double[::1] out_view = out
        cdef double R = self.R

        if A_view.shape[0] == 1:
            for i in range(n):
                out_view[i] = haversine(A_view[0, 0], A_view[0, 1], B_view[i, 0], B_view[i, 1], R)
        else:
            for i in range(n):
                out_view[i] = haversine(A_view[i, 0], A_view[i, 1], B_view[i, 0], B_view[i, 1], R)

        return out



