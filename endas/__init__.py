"""
Ensemble Data Assimilation package for Python.

The top-level :mod:`endas` package includes definitions of interfaces used by the library.
The :mod:`endas.algorithms` sub-package contains data assimilation algorithms, the
:mod:`endas.localization` sub-package contains classes for setting up localization of the analysis
update.

"""

from . import api
from .api import *

from . import cs
from . import cov
from . import algorithms
from . import arraycache


__all__ =  []
__all__ += api.__all__


# Patch the module names of members from imported sub-modules for Sphinx
for _xname in __all__:
    getattr(api, _xname).__module__ = __name__




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
