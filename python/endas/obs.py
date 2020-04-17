"""
Common observation operator implementations.

Currently only the matrix wrapper operator is implemented.
"""


import numpy as np

from . import ObservationOperator


class MatrixObservationOp(ObservationOperator):
    """
    Wraps (sparse) matrices as observation operators.

    The operator can be initialized with any NumPy/SciPy array or matrix, including sparse matrices. In fact, only
    the ``shape``, ``dot()`` and ``transpose()`` properties and methods are required.

    Args:
        H : Observation operator given as matrix.

    """

    def __init__(self, H):
        self._h = H


    @property
    def is_linear(self): return True


    @property
    def shape(self): return self._h.shape


    def localize(self, selected):
        return MatrixObservationOp(self._h[selected,:])


    def dot(self, x, out=None):
        if isinstance(self._h, np.ndarray): return self._h.dot(x, out=out)
        else:                               return self._h.dot(x)


    def to_matrix(self, force_dense=False, out=None):
        return self._h if not force_dense else np.asarray(self._h)



