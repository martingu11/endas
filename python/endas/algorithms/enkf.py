"""
Ensemble Kalman Filter and Smoother
===================================

The EnsembleKalmanFilter class implements Kalman Filter and Smoother based on an ensemble approximation of the error
covariance matrix.

EnsembleKalmanFilter relies on composition to implement different variants of the Kalman Filter/Smoother. Therefore,
its main purpose is to handle common functionality and domain (i.e. tile) -based localization of the analysis. The
actual Kalman update is then handled by subclasses of the `EnKFVariant` abstract base. This decoupling of execution and
the actual filter updates makes it easier to implement alternative approaches to how the sequential updates are
executed and in particular how localization is done.
"""

__all__ = [ 'EnsembleKalmanFilter', 'EnKFVariant', 'EnKF' ]

import numpy as np
from scipy import linalg

from python.endas import arraycache, ensemble
from python.endas.localization import DomainLocalization


from abc import ABCMeta, abstractmethod


class EnsembleKalmanFilter:
    """
    Ensemble Kalman Filter.

    This is a generic Ensemble Kalman Filter (EnKF) implementation that covers many of the EnKF variants. Ensemble
    Kalman Smoother (EnKS) is also implemented via the smoothing API. The filter/smoother allows both for global and
    localized analysis to be computed, the latter is always performed by partitioning the state space into a set of
    disjoint domains that are processed independently.

    Args:
        variant : Instance of the EnKFVariant class implementing the actual EnkF filter update
        ensemble_size : Size of the ensemble used
        loc_strategy : Localization strategy instance or ``None``.
        cov_inflation : Covariance inflation factor as a number greater or equal to 1.0. The default is 1.0,
                        meaning no covariance inflation is applied.
        lag : The lag length of the fixed-lag Kalman Smoother as a number of time steps. Zero Lag disables smoothing
        forgetting_factor: Forgetting factor applied by the smoother, the value must be less or equal to 1.0.
        cache : Large array cache instance or ``None``. If ``None`` is given, an internal memory cache will be used.

    See the :ref:`da-algorithms` page for information about implemented EnKF variants.
    """

    def __init__(self, variant, ensemble_size, loc_strategy=None,
                 cov_inflation=1.0, lag=0, forgetting_factor=1.0,
                 cache=None):

        self.localization_strategy = loc_strategy  # Localization strategy
        self.cov_inflation = cov_inflation  # Covariance inflation factor
        self.forgetting_factor = forgetting_factor

        assert isinstance(variant, EnKFVariant)
        self.variant = variant

        self._ensemblesize = ensemble_size
        self._lag = lag
        self._cache = arraycache.ArrayCache() if cache is None else cache
        self._smoother_data = []
        self._variant_initialized = False

        self.localize(loc_strategy)



    def localize(self, loc_strategy):
        """
        Applies localization to the Ensemble Kalman Filter.

        Args:
            loc_strategy: Localization strategy instance or ``None``.

        Returns:
            Nothing

        Calling ``localize()`` is equivalent to passing ``loc_strategy`` to the ``__init__`` method. Please note
        that the method must be called before ``smoother_begin()`` or ``begin_analysis()`` is called.
        """
        self._num_domains = 0
        self._loc_strategy = loc_strategy
        self._loc_statesize_sum = 0
        self._loc_statesize_limits = None

        if loc_strategy is None: return

        if isinstance(loc_strategy, DomainLocalization):
            self._num_domains = loc_strategy.ssp.num_domains
        else:
            raise TypeError("Unsupported localization strategy instance.")

        if self._num_domains == 0: return

        # We have local domains. We will store all local state vectors in a single 1-dimensional array and
        # pre-compute the array limits for each local state vector. In the future, we could consider using
        # multiple arrays if one becomes too large. Note: We do not rely on the sum of local state vector
        # lengths to equal the global state since there may be some padding applied to the local domains.

        self._loc_statesize_limits = np.empty((self._num_domains, 2), dtype=np.uint32)
        ssp = loc_strategy.ssp
        for di in range(self._num_domains):
            n = ssp.get_local_state_size(di)
            assert n > 0
            self._loc_statesize_limits[di, 0] = self._loc_statesize_sum
            self._loc_statesize_limits[di, 1] = n
            self._loc_statesize_sum+= n


    def forecast(self, model, A, Q, dt):
        n, N = A.shape
        if Q is not None: assert Q.shape == (n, n)

        # Move ensemble members forward in time and optionally perturb state vectors with model error
        model(A, dt)

        if Q is not None:
            QX = Q.random_multivariate_normal(N)
            QX = ensemble.center(QX, out=QX)
            np.add(A, QX, out=A)
        return A


    def smoother_begin(self, A0, t0):
        """
        Must be called once at the beginning when smoothing solution is being calculated.

        When the smoothing solution is calculated (i.e. with ``lag > 0``), this method must be called to initialize
        the smoother. This method does nothing if only the filtering solution is calculated (i.e. ``lag==0``).

        Args:
            A0 : Array representing the initial ensemble.
            t0 : Time corresponding to the initial time step.

        Returns:
            Nothing
        """
        if self._lag == 0: return

        self._X5 = None
        self._haveX5 = None

        self._Aa_enkf = None
        self._smoother_data = []
        self._cache.clear()

        n, N = A0.shape  # State and ensemble size

        if self._num_domains > 0:
            A0 = self._partition_state(A0)

        A0_handle = self._cache.put(A0)
        self._smoother_data.append((A0_handle, t0))


    def begin_analysis(self, A, t=None):
        """
        Initiates the analysis update step.

        At each time step, the analysis update must be initiated by a call to ``begin_analysis()``, followed by one or
        more calls to ``assimilate()``. After all observations have been assimilated for the current time step,
        ``end_analysis()`` must be called.

        Args:
            A : Array representing the forecast ensemble at current time step
            t : Time corresponding to the current filter/smoother step

        Returns:
            Nothing
        """

        n, N = A.shape  # State and ensemble size

        if not self._variant_initialized:
            self._variant_initialized = True
            self.variant.begin(n, N)

        self._n = n
        self._N = N
        self._t = t
        self._X5 = None
        self._haveX5 = None
        self._Af = A

        if self.cov_inflation != 1.0:
           ensemble.inflate(self._Af, self.cov_inflation, out=self._Af)

        if self._num_domains == 0:
            self._Aa_enkf = A
        else:
            self._Aa_enkf = self._partition_state(A)


    def assimilate(self, z, z_coords, H, R, executor=None):
        """
        Assimilates observations as part of the analysis update.

        Args:
            z        : Array of observations (values).
            z_coords : Observation coordinates. See below for more information.
            H        : Observation operator, must be an instance of :class:`endas.ObservationOperator`.
            R        : Observation error covariance, must be an instance of :class:`endas.CovarianceOperator`.

        Returns:
            Nothing

        The ``z_coords`` parameter contains "locations" of observations in the observation vector ``z`` and is only
        required when performing domain-based local analysis via ``localize()``. What can be passed for ``z_coords``
        depends on the state space partitioning (:class:`StateSpacePartitioning`) used for the localization, please
        check the documentation of the state space partitioning you are using for information. If the analysis is not
        localized, ``None`` can be passed for ``z_coords``.
        """

        n, N = self._n, self._N  # State and ensemble size
        m = len(z) if z is not None else 0  # Number of observations

        # Global analysis
        if self._num_domains == 0:
            if m > 0:
                Ag_data = self.variant.process_global_ensemble(self._Af, H)

                X5, X5s = self.variant.ensemble_transform(
                    self._Af, z, H, R,
                    Ag_data,
                    self.cov_inflation,
                    self.localization_strategy)
                self._Af = self._Af.dot(X5)
            else:
                X5s = None

            if self._X5 is not None: np.dot(self._X5, X5s, out=self._X5)
            else: self._X5 = X5s
            self._haveX5 = True
            self._Aa_enkf = self._Af


        # Localized analysis
        else:
            ls = self._loc_strategy
            ssp = self._loc_strategy.ssp

            if self._X5 is None:
                #Todo: Here we pack all X5 instances in a large array, which is wasteful if only a few domains
                #      are actually being updated with observations
                self._X5 = np.empty((self._num_domains, N, N))
                self._haveX5 = np.zeros(self._num_domains, dtype=np.bool)


            if m > 0:

                for di in range(self._num_domains):

                    # Collect observations for this domain...
                    local_zindexes, local_zdist = ssp.get_local_observations(di, z_coords, self._loc_strategy.taper_fn)
                    local_m = len(local_zindexes) if local_zindexes is not None else 0

                    # ... and assimilate if we have any
                    if local_m > 0:

                        # Local state vector
                        local_start, local_n = self._loc_statesize_limits[di]
                        local_A = self._Aa_enkf[local_start:local_start+local_n, :]
                        assert local_A.shape == (local_n, N)

                        # Construct localized versions of the observation operator and observation error covariance
                        local_H = ls.get_local_H(H, local_zindexes)
                        local_R = ls.get_local_R(R, local_zindexes, local_zdist)
                        local_z = z[local_zindexes]

                        # Before computing the ensemble transform, let the variant pre-calculate any data from the
                        # global ensemble it needs.
                        Ag_data = self.variant.process_global_ensemble(self._Af, local_H)

                        if executor is not None:
                            raise NotImplementedError()
                        else:
                            local_X5, local_X5s = self.variant.ensemble_transform(
                                local_A, local_z, local_H, local_R,
                                Ag_data,
                                self.cov_inflation,
                                self.localization_strategy)

                            if local_X5s is None: local_X5s = local_X5
                            assert local_X5.shape == (N, N)
                            assert local_X5s.shape == (N, N)

                            # Todo: Do this in-place?
                            local_Aa = local_A.dot(local_X5)
                            self._Aa_enkf[local_start:local_start + local_n, :] = local_Aa

                            if not self._haveX5[di]:
                                self._X5[di, :, :] = local_X5s
                                self._haveX5[di] = True
                            else:
                                self._X5[di, :, :] = np.dot(self._X5[di, :, :], local_X5s)


            # After we're done with the analysis, reconstruct the global ensemble as we are going to
            # need it either at next assimmilate() call or in end_analysis().
            for di in range(self._num_domains):
                local_start, local_n = self._loc_statesize_limits[di]
                ssp.put_local_state(di, self._Aa_enkf[local_start:local_start + local_n, :], self._Af)


    def end_analysis(self, on_smoother_result=None, result_args=tuple()):
        n, N = self._n, self._N  # State and ensemble size
        k = len(self._smoother_data)

        # First the lagged smoother by updating all previous EnKF states with the transformation matrix from the
        # current analysis update. For lag=0 this does nothing
        for j in range(k-1, k-1-self._lag, -1):
            if j < 0: break
            Aj_handle, tj = self._smoother_data[j]

            #print("Smoother up", k - j, Aj_handle, tj)

            Aj = self._cache.get(Aj_handle)
            assert Aj is not None

            Aj_is_result = j == k-self._lag  # Do we have the ready smoothing result?

            # For global analysis both A_enkf and X5 are the global ensemble and transform
            # arrays, respectively, so we only need to compute Aj*X5
            As = None
            eyeN = np.eye(N)

            if self._num_domains == 0:
                if self._X5 is not None:
                    if self.forgetting_factor != 1.0:
                        X4 = self._X5 - eyeN
                        X4 *= self.forgetting_factor
                        self._X5 = eyeN + X4
                    Aj = Aj.dot(self._X5, out=Aj)
                if Aj_is_result: As = Aj

            # For local analyses, Aj is an SxN array with S being the sum of state vector sizes for all local domains
            # (and it may be fo that S != the global state size `n` if there is some padding applied to domains!) and
            # N is the ensemble size. Likewise, X5 is an DxNxN array.
            else:
                As = np.empty((n, N)) if Aj_is_result else None

                for di in range(self._num_domains):

                    local_start, local_n = self._loc_statesize_limits[di]
                    local_Aj = Aj[local_start:local_start + local_n, :]
                    assert local_Aj.shape == (local_n, N)

                    if self._haveX5 is not None and self._haveX5[di]: # and self._last_analysis_j[d_i] <= j:
                        X5d = self._X5[di, :, :]
                        if self.forgetting_factor != 1.0:
                            X4 = X5d - eyeN
                            X4 *= self.forgetting_factor
                            X5d = eyeN + X4
                            self._X5[di, :, :] = X5d

                        local_Aj[:] = local_Aj.dot(X5d)

                    if Aj_is_result:
                        self._loc_strategy.ssp.put_local_state(di, local_Aj, As)


            if As is not None:
                if on_smoother_result is not None:
                    on_smoother_result(ensemble.mean(As), As, tj, result_args)
                #print("Remove Aj")
                self._cache.remove(Aj_handle)


        # Done with smoothing, if there was any. Store analysis result for this update step and return it.
        if self._lag == 0:
            if on_smoother_result is not None:
                on_smoother_result(ensemble.mean(self._Af), self._Af, self._t, result_args)
        else:
            Aa_handle = self._cache.put(self._Aa_enkf)
            #print("Add Aj", k)
            self._smoother_data.append((Aa_handle, self._t))

        self._X5 = None
        self._haveX5 = None
        self._Aa_enkf = None

        # Return the filtering solution
        return self._Af


    def smoother_finish(self, on_smoother_result, result_args=tuple()):
        n, N = self._n, self._N  # State and ensemble size
        k = len(self._smoother_data)

        As_buffer = np.empty((n, N)) if self._num_domains > 0 else None

        # First the lagged smoother by updating all previous EnKF states with the transformation matrix from the
        # current analysis update. For lag=0 this does nothing
        for j in range(k - 1, k - 1 - self._lag, -1):
            if j < 0: break
            Aj_handle, tj = self._smoother_data[j]

            # print("Smoother up", k - j, Aj_handle, tj)
            Aj = self._cache.get(Aj_handle)
            assert Aj is not None

            # For global analysis both A_enkf and X5 are the global ensemble and transform
            # arrays, respectively, so we only need to compute Aj*X5
            As = None
            if self._num_domains == 0:
                As = Aj
            else:
                for di in range(self._num_domains):
                    local_start, local_n = self._loc_statesize_limits[di]
                    local_Aj = Aj[local_start:local_start + local_n, :]
                    assert local_Aj.shape == (local_n, N)
                    self.localization_strategy.ssp.put_local_state(di, local_Aj, As_buffer)
                As = As_buffer

            if on_smoother_result is not None:
                on_smoother_result(ensemble.mean(As), As, tj, result_args)
                # print("Remove Aj")
            self._cache.remove(Aj_handle)


    def _partition_state(self, A, out=None):
        """
        Partitions the global state/ensemble into local state vectors. The local state vectors are all stored in a
        single 1-dimensional array.
        """
        n, N = A.shape
        if out is None:
            out = np.empty((self._loc_statesize_sum, N))
        else:
            assert out.shape == (self._loc_statesize_sum, N)

        ssp = self._loc_strategy.ssp
        for di in range(self._num_domains):
            d_start, d_n = self._loc_statesize_limits[di]
            out[d_start:d_start + d_n, :] = ssp.get_local_state(di, A)

        return out



class EnKFVariant(metaclass=ABCMeta):
    """
    Base class for various Ensemble Kalman Filter/Smoother variants.

    EnKFVariant implementations do most of the heavy-lifting of the EnsembleKalmanFilter
    by computing the transform applied to the forecast ensemble applied during the analysis
    update.
    """

    def begin(self, n, N):
        """
        Called once at the beginning of data assimilation.

        This method is called when :meth:`EnsembleKalmanFilter.begin()` method is called and implementations may place
        any initialization code here.
        Args:
            n: The size of the state vector
            N: The ensemble size
        """
        pass


    def process_global_ensemble(self, Ag, H):
        """
        Called before ``ensemble_transform`` to compute additional data from the global ensemble.

        Args:
            Ag : Array of shape (n, N) containing the global ensemble.
            H  : Observation operator instance. If analysis is localized, this is the local observation operator.


        This allows the implementation to calculate additional data from the global ensemble. The return value of this
        function is passed to ``ensemble_transform`` as the ``Ag_data`` argument. The default implementation returns
        ``None``.

        .. note: The purpose of the method is to avoid passing the global ensemble to the ``ensemble_transform`` method
                 when domain localization is used. In most cases, the (local) observation operator is applied to the
                 global ensemble, resulting in much smaller transformed ensemble. This is then passed to the
                 ``ensemble_transform`` method, which may be executed in a separate process or on a different machine
                 altogether.
        """
        return None


    @abstractmethod
    def ensemble_transform(self, A, z, H, R, Ag_data, inflation, lstrategy, out=None):
        """
        Computes analysis step as a transform applied to the ensemble.

        Assuming the analysis update is computed as A_a = A_f * T, this gives the
        transform `T`. The result must be an NxN matrix, where N is the ensemble
        size. If `out` is not None, the returned transform MUST be stored in the
        `out` instance in addition to being returned!

        Returns: NxN array where N is the ensemble size.
        """
        raise NotImplementedError




class EnKF(EnKFVariant):
    """
    Classic (stochastic) Ensemble Kalman Filter with perturbed observations.
    """

    def process_global_ensemble(self, Ag, H):
        AX = ensemble.to_anomaly(Ag)
        HAX = H.dot(AX)
        HA = H.dot(Ag)
        return HAX, HA


    def ensemble_transform(self, A, z, H, R, Ag_data, inflation, lstrategy, out=None):
        n, N = A.shape  # State and ensemble size
        m = H.shape[0]  # Number of observations

        HAX, HA = Ag_data

        # We have many observations or the observation error covariance operator
        # does not support the addition operator -> perform inversion using the
        # ensemble approximation to R
        if False:  # R.mc_only or m > N:
            K = None
        # R can be included explicitly in the inversion
        else:
            R_as_matrix = R.to_matrix(force_dense=False)
            HPHtR = HAX.dot(HAX.T)
            HPHtR+= (N - 1) * R_as_matrix

            K = linalg.solve(HPHtR, HAX, overwrite_a=True, overwrite_b=True).T
            #K = linalg.lstsq(HPHtR, HAX, overwrite_a=True, overwrite_b=True, cond=0.01)[0].T


        D = R.random_multivariate_normal(N)
        D = ensemble.center(D, out=D)

        D = np.add(D, z.reshape(-1, 1), out=D)
        D = np.subtract(D, HA, out=D)
        X4 = K.dot(D, out=out)

        X5 = np.eye(N) + X4
        #np.fill_diagonal(X4, X4.diagonal() + 1.0)  # X4 + I in-place
        return X5, None




