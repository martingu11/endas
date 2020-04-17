"""
Domain-based localization of the analysis update.
"""

__all__ = [
    'StateSpacePartitioning', 'GenericStateSpace1d', 'DomainLocalization'
]

from abc import ABCMeta, abstractmethod
import math
import numpy as np
import python.endas.localization



class StateSpacePartitioning(metaclass=ABCMeta):
    """
    Abstract base class for state space partitioning schemes.

    Implementations of this abstract base define how to partition the state space into local domains for localized
    analysis.
    """

    @property
    @abstractmethod
    def num_domains(self):
        """
        Returns the number of local domains in the state space partitioning.
        """
        pass


    @abstractmethod
    def get_local_observations(self, domain_id, z_coords, taper_fn):
        """
        Locates observations to be used for local analysis of the given domain.

        Args:
            domain_id    : Index of the domain for which observations should be retrieved. Domain indexes start at 0.
            z_coords     : Abstract description of the locations of observations in the observation vector.
            taper_fn     : Tapering function that defines the localization radius.
        Returns:
            Tuple ``(z_local, d)`` where ``z`` is a flat array of indices into the observation vector and ``d`` are
            the distances of each selected observation from the domain.

        The interpretation of ``z_coords`` is entirely up to the implementing class to decide and ``TypeError`` should
        be raised if an unsupported type is passed. Please see the documentation of implementing classes for information
        on what types of objects can be passed via ``z_coords``.
        """
        pass


    @abstractmethod
    def get_local_state_size(self, domain_id):
        """
        Returns the state vector size for the given domain.
        Args:
            domain_id : Index of the domain whose state size should be returned. Domain indexes start at 0.
        """
        pass

    @abstractmethod
    def get_local_state(self, domain_id, xg):
        """
        Returns local state vector (or ensemble of state vectors) for the given domain.

        Args:
            domain_id : Index of the domain whose local state vector be returned. Domain indexes start at 0.
            xg        : Global state vector or ensemble (ensemble members stored in columns).
        """
        pass

    @abstractmethod
    def put_local_state(self, domain_id, xl, xg):
        """
        Writes local state vector (or ensemble of state vectors) for the given domain back to the global state vector.

        Args:
            domain_id : Index of the domain whose local state vector be stored. Domain indexes start at 0.
            xl        : Local state vector or ensemble
            xg        : Global state vector or ensemble to copy the local data to
        """
        pass


class GenericStateSpace1d(StateSpacePartitioning):
    """
    Analysis localization assuming generic state space.

    This a generic form of domain-based analysis localization where each state space variable is assigned its own
    local domain. The local analysis update is therefore performed for each state variable individually. Under this
    scheme, the indexes of state vector elements also act as their respective "coordinates". Furthermore, observations
    are assumed to be directly related to state vector elements and their coordinates are also expressed as the indexes
    of the corresponding state vector elements.

    Args:
        n : The number of elements in the state vector

    .. note::
       Please note that this localization scheme is rather arbitrary and is unlikely to be suitable for any real
       geophysical system. This is especially true if non-trivial observation operators are used or in case of
       multivariate data assimilation. The primarily purpose of ``GenericStateSpace1d`` is for testing analysis
       localization on synthetic problems.

       The localization scheme is not suitable for large state spaces as the number of local domains is equal to the
       size of the state vector. This will make the localized analysis very slow for large *n*.
    """

    def __init__(self, n):
        assert n > 0
        self._n = n


    @property
    def num_domains(self):
        # Under this scheme each state variable has its own local domain
        return self._n

    def get_local_observations(self, domain_id, z_coords, taper_fn):
        assert isinstance(z_coords, np.ndarray)
        assert domain_id >= 0 and domain_id <= self.num_domains

        r = int(math.ceil(taper_fn.support_range))

        d_min = max(0, domain_id - r)
        d_max = min(self._n, domain_id + r)
        selected = np.where((z_coords > d_min) & (z_coords < d_max))[0]

        selcoords = z_coords[selected]
        dist = np.abs(np.subtract(selcoords, domain_id), dtype=np.double)
        return selected, dist


    def get_local_state_size(self, domain_id):
        assert domain_id >= 0 and domain_id <= self.num_domains
        return 1

    def get_local_state(self, domain_id, xg):
        assert domain_id >= 0 and domain_id <= self.num_domains
        return xg[domain_id]

    def put_local_state(self, domain_id, xl, xg):
        assert domain_id >= 0 and domain_id <= self.num_domains
        xg[domain_id] = xl





class DomainLocalization:
    """
    Domain-based localization of the analysis update.

    This is a general-purpose implementation of the :ref:`domain localization scheme<localization_domain>`.
    Domain localization scheme works by partitioning the global state space into a set of disjoint *local domains*.
    The analysis update is then done for each local domain individually and the local domains are assembled back
    into a global state/ensemble afterwards.

    Args:
        cs         : Coordinate system of the grid. Must be an instance of :class:`endas.localization.CoordinateSystem`
        taper_fn   : Taper function that defines the localization radius. Must be an instance of
                     :class:`endas.localization.TaperFn` or ``None`` (see notes below).
        on_local_R_needed : Called if this class does not know how to construct localized observation error covariance
                            operator from the global one.

    The tapering function (``taper_fn``) is required since it defines the localization radius. It is allowed to pass
    ``None`` during construction but the taper function must be set (via ``set_taper_fn()``) before the first analysis
    update at the latest. It is allowed to replace the taper function before each analysis update with a different one.

    One of the steps required from the domain localization is to construct localized versions of the global observation
    operator and the observation error covariance operator. The localized versions are subsets of the global operators
    that only represent observations used in the local analysis. Localization of the observation operator is done by the
    :meth:`DomainLocalization.get_local_H` method, which relies on ``ObservationOperator.localize()`` being implemented.
    """

    def __init__(self, ssp, taper_fn=None):
        assert ssp is not None
        assert isinstance(ssp, StateSpacePartitioning)
        self._ssp = ssp
        self._taper_fn = None
        if taper_fn is not None: self.set_taper_fn(taper_fn)


    @property
    def ssp(self):
        """
        The state space partitioning scheme.
        """
        return self._ssp


    @property
    def taper_fn(self):
        return self._taper_fn

    def set_taper_fn(self, taper_fn):
        """
        Sets the tapering function that defines the localization radius and how the influence of
        observations is truncated with distance.
        """
        assert taper_fn is not None
        assert isinstance(taper_fn, python.endas.localization.TaperFn)
        self._taper_fn = taper_fn


    def get_local_H(self, Hg, obs_used):
        """
        Returns localized observation operator for the given domain.

        Args:
            Hg       : Observation operator representing all observations
            obs_used : Flat array of indices into the observation vector pointing to observations to be used

        Returns:
            New ::class::`endas.ObservationOperator` instance or ``None`` if ``obs_used`` is an empty array.
        """
        if len(obs_used) == 0: return None
        Hl = Hg.localize(obs_used)
        assert Hl.shape == (len(obs_used), Hg.shape[1])
        return Hl


    def get_local_R(self, Rg, obs_used, obs_dist):
        """
        Returns localized observation error covariance operator for the given domain.

        Args:
            Rg       : Global observation error covariance operator
            obs_used : Flat array of indices into the observation vector pointing to selected observations
            obs_dist : Array of distances of selected observations from the local domain
            taper_fn : Tapering function instance
        """
        assert len(obs_used) == len(obs_dist)

        # We need to select a subset of the observation error covariance for the selected observations and apply the
        # tapering function to it.

        taper = np.ones(len(obs_dist), dtype=np.double)
        self._taper_fn.taper(taper, obs_dist, out=taper)

        return Rg.localize(obs_used, taper)


