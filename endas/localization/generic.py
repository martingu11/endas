"""
Domain localization on a generic state space.

"""

__all__ = ['LocalizationGeneric']

import math
import numpy as np
from endas.api import CoordinateSystem, DomainLocalization
from endas.localization.dl_base import DomainLocalizationBase

from endas import _get_cython_impl



class LocalizationGeneric(DomainLocalizationBase):
    """
    Implements domain localization of the analysis update assuming generic state space.

    This a generic form of domain -based analysis localization without any specific assumptions about the structure of
    the state space. Under this scheme each state variable belongs to its own local domain, i.e. the local analysis
    update is performed for each state variable individually.

    Args:
        n : The number of elements in the state vector
        cs : The coordinate system for measuring distances between state variables and observations. Must be an instance
             of :class:`endas.CoordinateSystem`

    .. note::
       Please note that this localization scheme is only suited to relatively small state spaces as the number of local
       domains is equal to the size of the state vector. This will make the localized analysis very slow for large *n*.

    """

    def __init__(self, n, cs):
        super().__init__(cs)
        assert n > 0
        self.n = n


    def generate_domains(self):
        # Generate domains. Since the domains are so trivial, there is nothing we need to store. It is enough to
        # have an index for each local domain and the rest we calculate as needed
        return [x for x in range(self.n)]

    def get_state(self, domain, xg):
        assert isinstance(domain, int)
        return xg[domain]

    def put_state(self, domain, xl, xg):
        assert isinstance(domain, int)
        xg[domain] = xl


