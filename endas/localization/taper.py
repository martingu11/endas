"""
Covariance tapering functions.

"""

from . import TaperFn
from endas import _get_cython_impl

import numpy as np

__all__ = ['GaspariCohn', 'Linear', 'Spherical']


class GaspariCohn(TaperFn):
    """
    Gaspari-Cohn covariance tapering function.

    The function is defined by

    .. math::

      G(r) =
      \\begin{cases}
        \\mbox{if } 0 \\leq r < 1 \\mbox{: } & 1 -      \\frac{5}{3}r^2 + \\frac{5}{8}r^3 + \\frac{1}{2}r^4 - \\frac{1}{4}r^5 \\\\
        \\mbox{if } 0 \\leq r < 2 \\mbox{: } & 4 - 5r + \\frac{5}{3}r^2 + \\frac{5}{8}r^3 - \\frac{1}{2}r^4 + \\frac{1}{12}r^5 - \\frac{2}{3r} \\\\
        \\mbox{otherwise : }                 & 0
      \end{cases}

    where :math:`r = \\vert i-j \\vert / L` and `i`, `j` are the (indexes of) two state variables,
    `L` is the distance scale factor (the localization radius). Therefore, the support range of the
    tapering function is `2L`.
    """

    def __init__(self, L):
        """
        Gaspari-Cohn covariance tapering function.

        Args:
            L : Correlation length of the tapering function. The tapering function will
                be zero at 2*L
        """
        assert L > 0
        self._L = L


    @property
    def support_range(self): return 2 * self._L

    def taper(self, x, d, out=None):
        assert isinstance(x, np.ndarray)
        assert isinstance(d, np.ndarray)
        assert x.ndim == 1
        assert d.ndim == 1
        assert len(x) == len(d)

        _cyimpl = _get_cython_impl()
        if _cyimpl is not None:
            if out is None: out = np.empty(len(x), dtype=np.double)
            return _cyimpl.taper_gc(x, d, out=out, L=self._L)
        # No Cython implementation available, fall-back on the NumPy version
        else:
            # Todo: Needs to be implemented!
            raise NotImplemented()



class Linear(TaperFn):
    """
    Linear covariance tapering function.
    The tapering function implements linear falloff towards zero at its support range `L`.
    """

    def __init__(self, L):
        assert L > 0
        self._L = L

    @property
    def supportrange(self): return self._L


    def taper(self, x, d, out=None):
        assert isinstance(x, np.ndarray)
        assert isinstance(d, np.ndarray)
        assert x.ndim == 1
        assert d.ndim == 1
        assert len(x) == len(d)

        _cyimpl = _get_cython_impl()
        if _cyimpl is not None:
            if out is None: out = np.empty(len(x), dtype=np.double)
            return _cyimpl.taper_linear(x, d, out=out, L=self._L)
        # No Cython implementation available, fall-back on the NumPy version
        else:
            # Todo: Needs to be implemented!
            raise NotImplemented()




class Spherical(TaperFn):
    """
    Spherical covariance tapering function.

    The taper function is defined by

    .. math::

      G(r) =
      \\begin{cases}
        \\left( 1 - \\left( \\frac{2}{3}\\frac{r}{L} - \\frac{1}{2}\\frac{r}{L}^3 \\right) \\right) & \\mbox{if } r < L \\\\
        0  & \\mbox{otherwise}
      \end{cases}

    where `r` is the distance (between two state variables) and `L` is the distance scale factor
    (the localization radius). Therefore, the support range of the tapering function is `L`.
    """

    def __init__(self, L):
        assert L > 0
        self._L = L

    @property
    def supportrange(self): return self._L


    def taper(self, x, d, out=None):
        assert isinstance(x, np.ndarray)
        assert isinstance(d, np.ndarray)
        assert x.ndim == 1
        assert d.ndim == 1
        assert len(x) == len(d)

        _cyimpl = _get_cython_impl()
        if _cyimpl is not None:
            if out is None: out = np.empty(len(x), dtype=np.double)
            return _cyimpl.taper_spherical(x, d, out=out, L=self._L)
        # No Cython implementation available, fall-back on the NumPy version
        else:
            # Todo: Needs to be implemented!
            raise NotImplemented()


