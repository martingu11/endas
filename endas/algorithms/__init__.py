"""
Data assimilation algorithms.

"""

from .kf import *
from .enkf_base import *
from .enkf import *
from .etkf import *


__all__ = []
__all__ += kf.__all__
__all__ += enkf_base.__all__
__all__ += enkf.__all__
__all__ += etkf.__all__


# Patch the module names of members from imported sub-modules for Sphinx
for _xname in __all__:
    globals()[_xname].__module__ = __name__



