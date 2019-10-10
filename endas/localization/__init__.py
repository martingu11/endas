"""
Localization of the analysis update.

"""

from abc import ABCMeta, abstractmethod

from . import domain
from .domain import *


__all__ = [
    'CoordinateSystem', 'TaperFn', 'SpatialQuery'
]
__all__.extend(domain.__all__)


#for x in domain.__all__: getattr(domain, x).__module__ = __name__


class TaperFn(metaclass=ABCMeta):
    """
    Covariance tapering function.

    Tapering functions are used for distance-based localization of the analysis update
    by adjusting the influence of observations based on their distance.
    """

    @property
    @abstractmethod
    def support_range(self):
        """
        Support range of the tapering function, i.e. the distance `d` at which the tapering
        coefficient becomes zero.
        """
        pass

    @abstractmethod
    def taper(self, x, d, out=None):
        """
        Tapers the vector `x` based on corresponding distances `d`.

        Tapering means multiplying each element `x_i` in the array `x` by weight `w(d_i)`,
        where `d_i` is the distance assigned to the element `x_i`.

        Args:
            x:   The array to taper
            d:   The distances assigned to elements in `x`
            out: Array where to store the result or None. If given, must have same shape as `x`

        Returns:
            The tapered array.
        """
        pass


class CoordinateSystem(metaclass=ABCMeta):
    """
    Coordinate system abstraction.

    Each coordinate system is characterized by its dimension (i.e. the number of components each
    coordinate takes) and must provide a distance metric via the `distance()` method.
    """

    @property
    @abstractmethod
    def ndim(self):
        """
        Returns the number of dimensions.
        """
        pass

    @property
    @abstractmethod
    def is_cartesian(self):
        """
        Returns ``True`` if this is a Cartesian coordinate system.
        """
        pass


    @abstractmethod
    def distance(self, A, B, out=None):
        """
        Implements distance metric for this coordinate system.

        The method computes distance between pairs of points in sets A and B. The dimension of each point must be equal
        to ``self.ndim`` and  are stored row-wise. Therefore, the sets A and B must be arrays of shape (n, ndim), where
        n is the number of points to process.

        Args:
            A, B: n x ndim arrays containing ``n`` coordinates in the sets ``A`` and ``B``
            out : Existing array of length ``n`` where the result should be stored, is possible. If ``None`` is given,
                  new array is allocated.

        Returns:
            Array of length `n` containing the computed distances.
        """
        pass



class SpatialQuery(metaclass=ABCMeta):
    """
    Lookup of observations or state vector elements by location.
    """

    @property
    def ndim(self):
        """
        Returns the number of spatial dimensions.
        This is a convenience shortcut to `self.cs.ndim`.
        """
        assert self.cs is not None
        assert isinstance(self.cs, CoordinateSystem)
        return self.cs.ndim



    @property
    @abstractmethod
    def cs(self):
        """Returns the coordinate system the query operates in."""
        pass


    @abstractmethod
    def range_query(self, lbound, ubound):
        """
        Looks up items based on a range of coordinates.

        The range is given by the lower and upper bounds, therefore all elements whose coordinates
        are in the interval [lbound_i, ubound_i], for all i in 1 .. `self.ndim`, are returned.

        Args:
            lbound: Array of length `self.ndim` of lower bound coordinates
            ubound: Array of length `self.ndim` of upper bound coordinates

        Returns:
            Indexes
        """
        pass



