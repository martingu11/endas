"""
Common covariance operator implementations.
"""

import math


import numpy as np
from numpy import random
from scipy import sparse, linalg

from . import CovarianceOperator

__all__ = [
    'to_matrix', 'DenseCovariance', 'DiagonalCovariance',
]


def to_matrix(C, force_dense=False, out=None):
    """
    Returns matrix representation of a covariance operator.

    If the passed instance is already a covariance matrix, it is returned as-is. Compatible matrix types are

    * ``numpy.ndarray``, ``numpy.matrix`` - returned as-is
    * subclass of ``scipy.sparse.spmatrix`` - returned as-is if ``force_dense=False``, otherwise converted to
      ``numpy.ndarray``.

    If ``C`` is an instance of :class:`endas.CovarianceOperator`, the result of ``C.to_matrix()`` is returned.
    Please note that not all covariance operators support conversion to explicit matrix form.

    Args:
        C : :class:`endas.CovarianceOperator` instance or NumPy array.
        force_dense : Value passed to :meth:`endas.CovarianceOperator.to_matrix()`.
        out : If specified, the array is used as the output buffer instead of allocating new array.
              The provided array must have the correct shape.

    Returns:
        Dense NumPy array (``ndarray``) or sparse SciPy array depending on ``C`` and ``force_dense``.

    Raises:
        TypeError: if C is of unsupported type
        ValueError: if C is a non-square matrix
        NotImplementedError: raised by :meth:`endas.CovarianceOperator.to_matrix()`
        MemoryError: if the covariance matrix array cannot be allocated

    """
    if isinstance(C, np.ndarray): ret = C
    elif isinstance(C, np.matrix): ret = C
    elif isinstance(C, sparse.spmatrix): ret = C.toarray(out=out)
    elif isinstance(C, CovarianceOperator): ret = C.to_matrix(force_dense, out=out)
    else:
        raise TypeError("C is not an instance of matrix or CovarianceOperator")

    if ret.shape[0] != ret.shape[1]:
        raise ValueError("C must be a square matrix")

    return ret



class DiagonalCovariance(CovarianceOperator):
    """
    Implements diagonal covariance matrix.

    The covariance operator is internally represented by a sparse matrix. Currently only the main diagonal
    is supported but other diagonals could be added, if needed. The operator supports all methods of
    :class:`CovarianceOperator`. The covariance can be instantiated with either the array of diagonal
    elements or the reciprocal (inverse) array. This can prevent numerical issues in situations where the inverse
    coefficients are near zero (thus leading to very large coefficients on the original diagonal) and if only
    the inverse coefficients are needed (such as when only ``solve()`` is called).

    Args:
          diag (nx1 array) : Vector of diagonal elements of the covariance matrix. Defines covariance matrix of
                             shape (n, n)

          invdiag (nx1 array) : Vector of inverse diagonal elements of the covariance matrix. Defines covariance
                                matrix of shape (n, n)

    Only one of the ``diag`` and ``invdiag`` arrays can be given.
    """

    def __init__(self, diag=None, invdiag=None):
        if diag is not None and invdiag is not None:
            raise ValueError("Only one of 'diag' and 'invdiag' arrays can be given.")
        elif diag is None and invdiag is None:
            raise ValueError("One of 'diag' or 'invdiag' arrays is needed.")
        elif diag is not None:
            self._diag_is_original = True
            self._Q = sparse.diags(diag)
            self._Qinv = sparse.diags(np.reciprocal(diag))
            self._sddiag = np.sqrt(diag)
        else:
            self._diag_is_original = False
            self._Q = sparse.diags(np.reciprocal(invdiag))
            self._Qinv = sparse.diags(invdiag)
            self._sddiag = np.sqrt(self._Q)


    @property
    def shape(self): return self._Q.shape

    @property
    def is_diagonal(self): return True

    @property
    def mc_only(self): return False

    def diagonal(self):
        """
        Returns the matrix diagonal as an array.

        Returns:
            Array of length ``self.shape[0]`` or ``self.shape[1]``.
        """
        return self._Q.diagonal()

    def inv_diagonal(self):
        """
        Returns the inverse of the matrix diagonal as an array.

        Returns:
            Array of length ``self.shape[0]`` or ``self.shape[1]``.
        """
        return self._Qinv.diagonal()


    def random_multivariate_normal(self, N=1):
        assert N >= 1
        n = self.shape[0]

        if N == 1:
            return random.randn(n) * self._sddiag
        else:
            rv = random.randn(n * N).reshape(n, N)
            return rv * self._sddiag.reshape(n, 1)

    def solve(self, b, overwrite_b=False):
        return self._Qinv.dot(b)


    def to_matrix(self, force_dense=False, out=None):
        return self._Q if not force_dense else np.asarray(self._Q.todense())


    def add_to(self, x):
        # Todo: custom add-to-diagonal ufunc?
        np.fill_diagonal(x, x.diagonal() + self._Q.diagonal())


    def localize(self, selected, taper):
        sel_diag = self._Q.diagonal()[selected]
        sel_invdiag = self._Qinv.diagonal()[selected]

        if taper is not None:
            assert len(selected) == len(taper)
            sel_invdiag = sel_invdiag * taper
            return DiagonalCovariance(invdiag=sel_invdiag)
        elif self._diag_is_original:
            return DiagonalCovariance(diag=sel_diag)
        else:
            return DiagonalCovariance(invdiag=sel_invdiag)




class DenseCovariance(CovarianceOperator):
    """
    Dense covariance matrix as CovarianceOperator.

    This is a trivial wrapper around a dense covariance matrix stored in a NumPy array. Use only on small problems
    when other representations don't fit.

    Args:
        C : Square NumPy matrix or array representing the covariance matrix.

    """

    def __init__(self, C):
        self._C = np.asarray(C)

    @property
    def shape(self): return self._C.shape

    @property
    def is_diagonal(self): return False

    @property
    def mc_only(self): return False


    def random_multivariate_normal(self, N=1):
        assert N >= 1
        n = self.shape[0]
        if N == 1:
            return random.multivariate_normal(np.zeros(n), self._C)
        else:
            return random.multivariate_normal(np.zeros(n), self._C, size=N).T


    def solve(self, b, overwrite_b=False):
        return linalg.solve(self._C, b, sym_pos=True, overwrite_b=overwrite_b)


    def to_matrix(self, **_ignored):
        return self._C

    def add_to(self, x):
        np.add(x, self._C, out=x)
