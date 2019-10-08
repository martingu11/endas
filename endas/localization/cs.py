"""
Coordinate Systems.
"""

__all__ = ['EuclideanCS']

import numpy as np

from endas.api import CoordinateSystem
from endas import _get_cython_impl


class EuclideanCS(CoordinateSystem):
    """
    Implements N-dimensional Cartesian coordinate system using Euclidean distance metric.
    """

    def __init__(self, ndim):
        """
        Implements N-dimensional coordinate system in Euclidean space.

        Args:
            ndim : Dimensionality of the coordinate system, must be >= 1
        """
        assert ndim >= 1
        self._ndim = ndim

    @property
    def ndim(self): return self._ndim

    @property
    def is_cartesian(self): return True

    def distance(self, A, B, out=None):
        assert A.shape == B.shape
        assert A.shape[1] == self.ndim
        n = A.shape[0]
        if n == 0: return np.empty(0)

        cyimpl = _get_cython_impl()

        if self.ndim == 1:
            return np.subtract(A, B, out=out)
        elif self.ndim == 2 and cyimpl is not None:
            if out is None: out = np.empty(n, dtype=np.double)
            cyimpl.cs_euclid_distance_2d(A, B, out)
            return out
        elif self.ndim == 3 and cyimpl is not None:
            if out is None: out = np.empty(n, dtype=np.double)
            cyimpl.cs_euclid_distance_3d(A, B, out)
            return out
        # Cython variant not available, revert to NumPy which will require several passes
        else:
            d = np.subtract(A, B)
            d = np.square(d, out=d)
            d = np.sum(d, axis=1, out=out)
            d = np.sqrt(d, out=d)
            return d


class LatLonCS(CoordinateSystem):
    """
    Implements coordinate system on a perfect sphere with coordinates of any point expressed as latitude and longitude.

    The currently implemented coordinate system is strictly two-dimensional, i.e. there is no vertical component.
    The coordinates are assumed to be in degrees.

    Args:
        R : The radius of the great circle on which the distance is calculated. The default value of 6371 km corresponds
            to the mean Earth radius (R1) as defined by the International Union of Geodesy and Geophysics.

    This class is a simple implementation of a polar coordinate system assuming a perfect sphere of radius *R* and uses
    the `Haversine formula <https://en.wikipedia.org/wiki/Haversine_formula>`_ to compute the spherical distance.
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

        cyimpl = _get_cython_impl()

        if cyimpl is not None:
            if out is None: out = np.empty(n, dtype=np.double)
            cyimpl.cs_latlon_distance(A, B, out, self.R)
            return out
        # Cython variant not available, revert to NumPy which will require several passes
        else:
            raise NotImplementedError()



