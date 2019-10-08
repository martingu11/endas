"""
Abstract base for comain localization implementations.

"""

__all__ = ['DomainLocalizationBase']

import math
import numpy as np
from endas.api import CoordinateSystem, DomainLocalization
from endas import _get_cython_impl



class DomainLocalizationBase(DomainLocalization):
    """
    Abstract base class for domain localization implementations.

    This covers the common functionality that classes implementing domain localization need to deal with.
    """

    def __init__(self, cs):
        assert cs is not None
        assert isinstance(cs, CoordinateSystem)
        self.cs = cs

    @property
    def ndim(self): return self.cs.ndim


