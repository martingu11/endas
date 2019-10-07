"""
EnDAS
=====

Ensemble Data Assimilation package for Python.

"""

from . import api
from .api import *

from . import cs
from . import cov
from . import algorithms
from . import arraycache
from . import localization


__all__ =  []
__all__ += api.__all__


# Try loading the Cython implementation extension module. It may not have been compiled,
# which is fine, we have fall-back implementations in place.
_HAVE_CYTHON_IMPL = False
try:
    import endas._cython_impl as _cyimpl
    _HAVE_CYTHON_IMPL = True
except:
    pass

def _get_cython_impl():
    return _cyimpl if _HAVE_CYTHON_IMPL else None
