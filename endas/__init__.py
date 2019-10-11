"""
Ensemble Data Assimilation package for Python.

"""

from abc import ABCMeta, abstractmethod

__all__ = ['ObservationOperator', 'CovarianceOperator']


class ObservationOperator(metaclass=ABCMeta):
    """
    Abstract base class for observation operators.

    Observation operators provide mapping from the state space to the observation space. In other words, the operator
    transforms the state vector to the corresponding set of observations that we would expect.

    The operator may be implemented as a (typically sparse) matrix but more abstract implementations are possible too.
    """

    @property
    @abstractmethod
    def is_linear(self):
        """
        True if the operator is linear.
        """
        pass

    @property
    @abstractmethod
    def shape(self):
        """
        Returns the "shape" of the operator as a tuple ``(k, n)``, where ``k`` is the number of observations this
        operator represents and ``n`` is the state vector size.
        """
        pass


    @abstractmethod
    def localize(self, selected):
        """
        Returns new ``ObservationOperator`` for a subset of observations.

        Args:
            selected : Array of indexes of observations to select

        Returns:
            An instance of ``ObservationOperator`` of shape (s, n), where ``s = len(selected)``.

        Raises:
              NotImplementedError: if the operation is not supported by this operator

        """
        raise NotImplementedError


    @abstractmethod
    def dot(self, x, out=None):
        """
        Applies the observation operator to an ensemble or state vector.

        This implements the matrix product :math:`\mathbf{H}x` if ``H`` is linear, or more generally
        :math:`\mathscr{H}(x)` for non-linear observation operators. Here :math:`x` is an nxN array (i.e. either the
        ensemble or state vector).

        Args:
            x  : Vector or matrix of shape nxN to which the operator is applied to, where :math:`N \geq 1`.
            out : If supplied, it is used as a buffer for the output instead of allocating new array. The passed array
                  must be of shape kxN, where `k=self.shape[0]``.

        Returns:
            The result as kxN array, where ``k=self.shape[0]``.
        """
        pass


    @abstractmethod
    def to_matrix(self, force_dense=False, out=None):
        """
        Returns the matrix form of the operator, if available.

        Args:
           force_dense : If ``True``, the returned matrix is always a dense Numpy array. Otherwise a sparse matrix may
                         be returned, if suitable for the observation operator.
           out : If supplied, it is used as a buffer for the output instead of allocating new array. The passed array
                 must be same shape as `self.shape``.

        Raises:
              NotImplementedError: if the operation is not supported by this covariance operator

        Note:
            Please note that this should only be used on sparse observation operators and exists primarily so that the
            (Extended) Kalman Filter can be implemented.
        """
        raise NotImplementedError()



class CovarianceOperator(metaclass=ABCMeta):
    """
    Abstract representation of a covariance matrix.

    The base class defines interface required by covariance matrix implementations. See the ``endas.cov`` module for
    concrete implementations.

    """

    @property
    @abstractmethod
    def shape(self):
        """
        Returns the shape of the covariance matrix represented by this operator. Since this is a covariance matrix,
        the shape is always square
        """
        pass


    @property
    @abstractmethod
    def is_diagonal(self):
        """
        Returns ``True`` if the covariance is represented by a diagonal matrix.

        The default implementation returns ``False``.

        Note:
            If True, calling ``to_matrix(False)`` must return an instance of ``scipy.sparse.dia_matrix``, i.e. a
            sparse matrix in `diagonal` format.
        """
        return False



    @property
    @abstractmethod
    def mc_only(self):
        """
        Returns ``True`` if this covariance operator only supports Monte-Carlo sampling
        via the `random_multivariate_normal()`.

        The default implementation returns ``True`` as this is the minimal functionality a covariance operator must
        support. Implementations that return `False` will throw `NotImplementedError` from any method other than
        `random_multivariate_normal()`.
        """
        return True


    @abstractmethod
    def to_matrix(self, force_dense=False, out=None):
        """
        Returns the covariance matrix.

        Args:
           force_dense : If True, the returned matrix is always dense. Otherwise a sparse matrix may be returned,
                         if suitable for the covariance operator.
           out : If specified, the array is used as the output buffer instead of allocating new array. The provided
                 array must have the correct shape.

        The default implementation raises NotImplementedError.

        Raises:
            NotImplementedError: if the operation is not supported
            MemoryError: if the covariance matrix array cannot be allocated

        Note:
            Please note that this should only be used on special cases such as diagonal covariance matrices
            (see :func:`is_diagonal`) or for visualization and debugging. The call is likely to run out of memory
            when use on a large covariance operator.
        """
        raise NotImplementedError()



    @abstractmethod
    def random_multivariate_normal(self, N=1):
        """
        Implements generation of a random sample from multivariate Normal distribution with zero mean and covariance
        given by this CovarianceMatrix instance.

        Args:
          N : The number of independent samples to draw.

        Returns:
          nxN array where `n` is the state space size (i.e. `self.shape[0]`).
        """
        pass


    @abstractmethod
    def solve(self, b, overwrite_b=False):
        """
        Solves the system :math:`\mathbf{C}x = b`, where :math:`\mathbf{C}` is the covariance matrix represented
        by this ``CovarianceOperator`` instance.

        Args:
          b (n x m array): Data for the right hand side, where n is the size of the
                           covariance matrix and :math:`m \geq 1`

        Returns:
          The solution `x`. The shape of the array is equivalent to the shape of `b`.
        """
        raise NotImplementedError()


    @abstractmethod
    def add_to(self, x):
        """
        Sums this covariance matrix and `x` and stores the result in `x`. Please note that
        this may not be supported by all covariance matrix implementations.

        Generally, only covariance matrices that can afford to be represented explicitly (small
        matrices, diagonal matrices) implement this. In other cases ``NotImplementedError`` is
        raised.
        """
        raise NotImplementedError()


    def localize(self, selected, taper):
        """
        Returns new ``CovarianceOperator`` for a subset of the original operator.

        The default implementation raises ``NotImplementedError``.

        Args:
            selected : Array of indexes of observations to select

        Returns:
            An instance of ``ObservationOperator`` of shape (s, n), where ``s = len(selected)``.

        Raises:
              NotImplementedError: if the operation is not supported by this operator
        """
        raise NotImplementedError
