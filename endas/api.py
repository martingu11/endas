#
# This module contains all public interfaces of the EnDAS library. All public members of
# this module are going to be pulled into the top-level `endas` package namespace.
#

__all__ = [
    'CoordinateSystem', 'ObservationOperator', 'CovarianceOperator',
    'TaperFn', 'SpatialQuery', 'DomainLocalization'
    ]

from abc import ABCMeta, abstractmethod
import numpy as np


class ObservationOperator(metaclass=ABCMeta):
    """
    Abstract base class for observation operators.

    Observation operators provide mapping from the state space to the observation space.
    In other words, the operator transforms the state vector to the corresponding set of
    observations that we would expect.

    The operator may be implemented as a (typically sparse) matrix but more abstract
    implementations are possible too.
    """


    @property
    @abstractmethod
    def islinear(self):
        """
        True if the operator is linear.
        """
        pass

    @property
    @abstractmethod
    def shape(self):
        """
        Returns the "shape" of the operator as a tuple ``(k, n)``, where `k` is the number of
        observations and `n` is the state vector size.
        """
        pass

    @abstractmethod
    def slice(self, r):
        """
        Returns a "slice" of the operator, i.e. an identical operator restricted to the observations
        in `r`. The returned operator should therefore have shape ``(len(r), n)``.

        The returned value must be an instance of ``ObservationOperator``.
        """
        raise NotImplementedError

    @abstractmethod
    def slice_range(self, rmin, rmax):
        """
        Returns a "slice" of the operator, i.e. an identical operator restricted to the
        range of observations [rmin, rmax). The returned operator should therefore have
        shape ``(rmax-rmin, n)``.

        The returned value must be an instance of ``ObservationOperator``.
        """
        raise NotImplementedError



    @abstractmethod
    def dot(self, x, out=None):
        """
        Implements :math:`\mathbf{H}x` (or :math:`\mathscr{H}(x)`).

        Args:
            x  (n x a array) : Vector or matrix to which the operator is applied to, where :math:`a \geq 1`.

        Returns:
            The result as n x a array.
        """
        pass


    @abstractmethod
    def adjdot(self, x, out=None):
        """
        Implements adjoint of the operator, i.e. :math:`\mathbf{H}^{\mathsf{T}} x` (or
        :math:`\mathscr{H}^{\mathsf{T}}(x)`).

        Args:
            x  (k x a array) : Vector or matrix to which the operator is applied to, where :math:`a \geq 1`.

        Returns:
            The result as n x a array.
        """
        pass


    def tomatrix(self, forceDense=False):
        """
        Returns the matrix form of the operator.

        Args:
           forceDense : If True, the returned matrix is always dense. Otherwise a sparse
                        matrix may be returned, if suitable for the observation operator.

        The default implementation raises NotImplementedError.

        Note:
            Please note that this should only be used on sparse observation operators and
            primarily exists so that the Extended Kalman Filter can be implemented.
        """
        raise NotImplementedError()



class CovarianceOperator(metaclass=ABCMeta):
    """
    Abstract representation of a covariance matrix.

    The base class defines interface required by covariance matrix implementations.
    See the ``endas.cov`` module for concrete implementations.

    """

    @property
    @abstractmethod
    def shape(self):
        """
        Returns the shape of the covariance matrix represented by this operator. Since this is a
        covariance matrix, the shape is always square
        """
        pass


    @property
    def is_diagonal(self):
        """
        Returns ``True`` if the covariance is represented by a diagonal matrix.

        The default implementation returns ``False``.

        Note:
            If True, calling ``to_matrix(False)`` must return an instance of ``scipy.sparse.dia_matrix``,
            i.e. a sparse matrix in `diagonal` format.
        """
        return False



    @property
    def mc_only(self):
        """
        Returns ``True`` if this covariance operator only supports Monte-Carlo sampling
        via the `random_multivariate_normal()`.

        The default implementation returns ``True`` as this is the minimal functionality
        a covariance operator must support. Implementations that return `False` will throw
        `NotImplementedError` from any method other than `random_multivariate_normal()`.
        """
        return True



    def to_matrix(self, force_dense=False):
        """
        Returns the covariance matrix.

        Args:
           force_dense : If True, the returned matrix is always dense. Otherwise a sparse
                         matrix may be returned, if suitable for the covariance operator.

        The default implementation raises NotImplementedError.

        Note:
          Please note that this should only be used on special cases such as diagonal
          covariance matrices (see :func:`is_diagonal` or for visualization and debugging.
          when the matrix dimension is small enough.
        """
        raise NotImplementedError()

    @classmethod
    def as_matrix(cls, instance, forceDense=False):
        if isinstance(instance, np.ndarray): return instance
        elif isinstance(instance, np.matrix): return instance
        elif isinstance(instance, CovarianceOperator): return instance.tomatrix(forceDense)



    @abstractmethod
    def random_multivariate_normal(self, N=1):
        """
        Implements generation of a random sample from multivariate Normal distribution
        with zero mean and covariance given by this CovarianceMatrix instance.

        Args:
          N : The number of independent samples to draw.

        Returns:
          nxN array where `n` is the state space size (i.e. `self.shape[0]`).
        """
        pass


    @abstractmethod
    def solve(self, b, overwrite_b=False):
        """
        Solves the system :math:`\mathbf{C}x = b`, where :math:`\mathbf{C}`
        is the covariance matrix represented by this CovarianceMatrix instance.

        Args:
          b (n x m array): Data for the right hand side, where n is the size of the
                           covariance matrix and :math:`m \geq 1`

        Returns:
          The solution `x`. The shape of the array is equivalent to the shape of `b`.
        """
        raise NotImplementedError()


    def add_to(self, x):
        """
        Sums this covariance matrix and `x` and stores the result in `x`. Please note that
        this may not be supported by all covariance matrix implementations.

        Generally, only covariance matrices that can afford to be represented explicitly (small
        matrices, diagonal matrices) implement this. In other cases ``NotImplementedError`` is
        raised.
        """
        raise NotImplementedError()
    #
    # def sqrt_addto(self, x):
    #     """
    #     Sums thie square root of this covariance matrix and `x` and stores the result in `x`.
    #     Please note that this may not be supported by all covariance matrix implementations.
    #
    #     Generally, only covariance matrices that can afford to be represented explicitly (small
    #     matrices, diagonal matrices) implement this. In other cases ``NotImplementedError`` is
    #     raised.
    #     """
    #     raise NotImplementedError()



class SequentialEnsembleAlgorithm(metaclass=ABCMeta):
    """
    Abstract base class for sequential ensemble-based data assimilation algorithms.

    SequentialAlgorithm has API for both sequential filtering and smoothing though the latter may
    not be implemented by all algorithms. The filtering API is always implemented.

    """

    @abstractmethod
    def begin(self, E0, t0):
        """
        Args:
            E0:
            t0:

        Returns:
        """
        pass


class TaperFn(metaclass=ABCMeta):
    """
    Tapering function.

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
    coordinate takes) and must provide a distance metric via the `distance()` medhod.
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




class DomainLocalization(metaclass=ABCMeta):
    """
    Localization of the analysis update by partitioning of the state space to disjoint sub-domains
    for which local analysis is done.
    """

    @abstractmethod
    def set_obscoords(self, coords):
        """
        Sets the observation coordinates.


        This is a good place for implementations to do something with the coordinates for
        efficient retrieval of local observations, such as creating a spatial index.
        """
        pass


    @abstractmethod
    def generate_domains(self):
        """
        Returns sequence of distinct domains over the state space.

        What comprises a 'domain' is entirely to the implementing class. The returned object may
        be any iterable sequence, including generator objects or expressions. The contents of the
        sequence are also specific to the implementing class and are to be passed as they are into
        `get_observation()` etc.
        """
        pass


    @abstractmethod
    def get_obs_data(self, domain, Hg, Rg):
        """
        Returns local observation data for the given domain.

        Returns tuple (z, H, R)
        """
        pass


    def get_H(self, domain, Hg, obsused):
        """
        Returns localized observation operator for the given domain.

        Args:
          Hg : Global observation operator
        """
        pass


    def get_R(self, domain, Rg, obsused, dist, taperfn):
        """
        Returns localized observation error covariance operator for the given domain.

        Args:
            Rg : Global observation error covariance operator
        """
        pass


    @abstractmethod
    def get_state_size(self, domain):
        """
        Returns the state vector size for the given domain.
        Args:
            domain : Domain instance whose state size should be returned
        """
        pass

    @abstractmethod
    def get_state(self, domain, xg):
        """
        Returns local state vector (or ensemble of state vectors) for the given domain.

        Args:
          xg : Global state vector or ensemble (ensemble members stored in columns)
        """
        pass


    @abstractmethod
    def put_state(self, domain, xl, xg):
        """
        Writes local state vector (or ensemble of state vectors) for the given domain back to
        the global state vector.

        Args:
          xl : Local state vector or ensemble
          xg : Global state vector or ensemble to copy the local data to
        """
        pass





