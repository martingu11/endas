"""
Spatial localization of the analysis update.
"""

__all__ = [
    'StateSpacePartitioning', 'GenericStateSpace1d', 'DomainLocalization'
]

from abc import ABCMeta, abstractmethod
import math
import numpy as np
import endas.cov as cov
import endas.localization



class StateSpacePartitioning(metaclass=ABCMeta):
    """
    Abstract base class for state space partitioning schemes.

    Implementations of this abstract base define how to partition the state space into local domains for localized
    analysis.
    """

    @property
    @abstractmethod
    def ndim(self):
        """
        Number of dimensions over which the partitioning is done.
        """
        pass

    @abstractmethod
    def generate_domains(self):
        """
        Returns sequence of distinct domains over the state space.

        The returned object may be any iterable sequence, including generator objects or expressions,
        containing/yielding the individual local domains. What comprises a "local domain" is entirely to the
        implementing class. Therefore, the domains as returned by this call are to be passed as they are to the
        ``get_local_*`` and other methods of this class to identify the local domain in question.
        """
        pass


    @abstractmethod
    def get_local_observations(self, domain, z_coords, taper_fn, distances=True):
        """
        Locates observations to be used for local analysis of the given domain.

        Args:
            domain    : Identifies domain for which observations should be retrieved.
            z_coords  : Abstract description of the locations of observations in the observation vector.
            taper_fn  : Tapering function that defines the localization radius.
            distances : If ``True``, distances of the returned observations from the domain are also
                        returned.
        Returns:
            Tuple ``(z_local, d)`` where ``z`` is a flat array of indices into the observation vector and ``d`` are
            the distances of each selected observation from the domain (or ``None`` if ``distances=False`` was passed).

        The interpretation of ``z_coords`` is entirely up to the implementing class to decide and ``TypeError`` should
        be raised if an unsupported type is passed. Please see the documentation of implementing classes for information
        on what types of objects can be passed via ``z_coords``.
        """
        pass


    @abstractmethod
    def get_local_state_size(self, domain):
        """
        Returns the state vector size for the given domain.
        Args:
            domain : Domain instance whose state size should be returned
        """
        pass

    @abstractmethod
    def get_local_state(self, domain, xg):
        """
        Returns local state vector (or ensemble of state vectors) for the given domain.

        Args:
          xg : Global state vector or ensemble (ensemble members stored in columns)
        """
        pass

    @abstractmethod
    def put_local_state(self, domain, xl, xg):
        """
        Writes local state vector (or ensemble of state vectors) for the given domain back to
        the global state vector.

        Args:
          xl : Local state vector or ensemble
          xg : Global state vector or ensemble to copy the local data to
        """
        pass


class GenericStateSpace1d(StateSpacePartitioning):
    """
    Analysis localization assuming generic state space.

    This a generic form of domain-based analysis localization where each state space variable is assigned its own
    local domain. The local analysis update is therefore performed for each state variable individually. Under this
    scheme the indexes of state vector elements also act as their respective "coordinates". Furthermore, observations
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

    #@property
    #def ndim(self): return 1

    def generate_domains(self):
        # Generate domains. Since the domains are so trivial, there is nothing we need to store. It is enough to
        # have an index for each local domain and the rest we calculate as needed
        return [x for x in range(self._n)]


    def get_local_observations(self, domain, z_coords, taper_fn, distances=True):
        assert isinstance(z_coords, np.ndarray)
        r = int(math.ceil(taper_fn.support_range))
        d_min = max(0, domain - r)
        d_max = max(self._n, domain + r)

        selected_z = np.where((z_coords > d_min) & (z_coords <= d_max))[0]

        if distances:
            selected_d = np.subtract(selected_z, domain)
            selected_d = np.abs(selected_d, out=selected_d)
        else:
            selected_d = None

        return selected_z, selected_d

    def get_local_state_size(self, domain):
        return 1

    def get_local_state(self, domain, xg):
        assert isinstance(domain, int)
        return xg[domain]

    def put_local_state(self, domain, xl, xg):
        assert isinstance(domain, int)
        xg[domain] = xl





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


    def set_taper_fn(self, taper_fn):
        """
        Sets the tapering function that defines the localization radius and how the influence of
        observations is truncated with distance.
        """
        assert taper_fn is not None
        assert isinstance(taper_fn, endas.localization.TaperFn)
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


