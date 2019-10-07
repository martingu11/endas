"""
Coordinate Systems.
"""

__all__ = ['EuclideanCS']

import numpy as np
from .api import CoordinateSystem


class EuclideanCS(CoordinateSystem):
    """
    Implements N-dimensional coordinate system in Euclidean space.
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


    def distance(self, A, B, out=None):
        assert A.shape == B.shape
        assert A.shape[1] == self.ndim
        if A.shape[0] == 0: return np.empty(0)


        # Cython variant not available, revert to NumPy which will require several passes
        d = np.subtract(A, B)
        d = np.square(d,out=d)
        d = np.sum(d, axis=1, out=out)
        d = np.sqrt(d, out=d)
        return d



