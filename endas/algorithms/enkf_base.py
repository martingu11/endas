"""
Ensemble Kalman Filter and Smoother
===================================

The EnsembleKalmanFilter class implements Kalman Filter and Smoother based on an ensemble
approximation of the error covariance matrix.


EnsembleKalmanFilter relies on composition to implement different variants of the Kalman
Filter/Smoother. Therefore, it's main purpose is to handle common functionality and
domain (i.e. tile) -based localization of the analysis. The actual Kalman update is then
handled by subclasses of the `EnKFVariant` abstract base. This decoupling of execution and
the actual filter updates makes it easier to implement alternative approaches to how the
sequential updates are executed and in particular how localization is done.

"""

__all__ = [ 'EnsembleKalmanFilter', 'EnKFVariant' ]

import math
import numpy as np

from endas import ensemble, CovarianceOperator
from endas import DomainLocalization
from endas import arraycache


from abc import ABCMeta, abstractmethod


class EnKFVariant(metaclass=ABCMeta):
    """
    Base class for various Ensemble Kalman Filter/Smoother variants.

    EnKFVariant implementations do most of the heavy-lifting of the EnsembleKalmanFilter
    by computing the transform applied to the forecast ensemble applied during the analysis
    update.
    """

    @abstractmethod
    def compute_ensemble_transform(self, A_global, A, z, H, R, inflation, lstrategy, out=None):
        """
        Computes analysis step as a transform applied to the ensemble.

        Assuming the analysis update is computed as A_a = A_f * T, this gives the
        transform `T`. The result must be an NxN matrix, where N is the ensemble
        size. If `out` is not None, the returned transform MUST be stored in the
        `out` instance in addition to being returned!

        Args:
            A :
            AX :
            HPHtR_inv :
            t :
            z :
            H :
            R :
            out :

        Returns: NxN array where N is the ensemble size.
        """
        raise NotImplementedError



class EnsembleKalmanFilter:
    """
    Ensemble Kalman Filter.

    This is a generic Ensemble Kalman Filter (EnKF) implementation that covers many of the
    EnKF variants. Ensemble Kalman Smoother (EnKS) is also implemented via the smoothing
    API. The filter/smoother allows both for global and localized analysis to be computed,
    the latter is always performed by partitioning the state space into a set of disjoint
    domains that are processed independently.

    See the `EnKFVariant` and its subclasses for implemented EnKF/EnKS variants.
    """

    def __init__(self, variant, ensemblesize, loc_strategy=None, cov_inflation=1.0, lag=0,
                 forgetting_factor=1.0):
        """
        Constructor.
        """
        self.localization_strategy = loc_strategy  # Localization strategy
        self.cov_inflation = cov_inflation  # Covariance inflation factor
        self.forgetting_factor = forgetting_factor

        assert isinstance(variant, EnKFVariant)
        self.variant = variant

        self.do_checks = True
        self.tempdir = None

        self._ensemblesize = ensemblesize
        self._haveArrays = False

        self._lag = lag
        self._cache = arraycache.ArrayCache()
        self._smoother_data = []

        self.localize(loc_strategy)



    def localize(self, loc_strategy):
        self._domains = []
        if loc_strategy is not None and isinstance(loc_strategy, DomainLocalization):
            self._domains = loc_strategy.generate_domains()

        self._loc_statesize_sum = 0
        self._loc_statesize_limits = None

        num_domains = len(self._domains)
        if num_domains == 0: return

        self._loc_statesize_limits = np.empty((num_domains, 2), dtype=np.uint32)
        for di, d in enumerate(self._domains):
            n = loc_strategy.get_state_size(d) if d is not None else 0
            #print (d, n)
            assert d is None or n > 0

            self._loc_statesize_limits[di, 0] = self._loc_statesize_sum
            self._loc_statesize_limits[di, 1] = n
            self._loc_statesize_sum+= n


    def forecast(self, model, dt, A, Q):
        n, N = A.shape
        if Q is not None: assert Q.shape == (n, n)

        # Move ensemble members forward in time and optionally perturb state vectors with model error
        model(A, dt)

        if Q is not None:
            QX = Q.random_multivariate_normal(N)
            QX = ensemble.center(QX, out=QX)
            np.add(A, QX, out=A)
        return A


    def begin(self, A0, t0):
        self._X5 = None
        self._haveX5 = None

        self._Aa_enkf = None
        self._smoother_data = []
        self._cache.clear()

        n, N = A0.shape  # State and ensemble size
        num_domains = len(self._domains)

        if num_domains > 0:
            A0_decomposed = np.empty((self._loc_statesize_sum, N))
            for d_i, d in enumerate(self._domains):
                if d is None: continue
                d_start = self._loc_statesize_limits[d_i][0]
                d_n = self._loc_statesize_limits[d_i][1]
                A0_decomposed[d_start:d_start + d_n, :] = self.localization_strategy.get_state(d, A0)

            A0 = A0_decomposed

        A0_handle = self._cache.put(A0)
        self._smoother_data.append((A0_handle, t0))

        self._last_analysis_j = np.zeros(num_domains, dtype=np.uint) if num_domains > 0 else 0




    def begin_analysis(self, A, t=None):
        n, N = A.shape  # State and ensemble size
        num_domains = len(self._domains)

        self._n = n
        self._N = N
        self._t = t
        self._X5 = None
        self._haveX5 = None
        self._Af = A

        # The number of assimilated observations (in each domain) during this analysis update
        self._num_z = np.zeros(num_domains, dtype=np.uint) if num_domains > 0 else 0

        #if self.cov_inflation != 1.0:
        #   ensemble.inflate(self._Af, self.cov_inflation, out=self._Af)

        if num_domains == 0:
            self._Aa_enkf = A
        else:
            self._Aa_enkf = np.empty((self._loc_statesize_sum, N))
            for d_i, d in enumerate(self._domains):
                if d is None: continue
                d_start = self._loc_statesize_limits[d_i][0]
                d_n = self._loc_statesize_limits[d_i][1]
                self._Aa_enkf[d_start:d_start + d_n, :] = self.localization_strategy.get_state(d, A)


    def assimilate(self, z, H, R):
        n, N = self._n, self._N  # State and ensemble size
        m = len(z) if z is not None else 0  # Number of observations
        num_domains = len(self._domains)

        # Global analysis
        if num_domains == 0:
            if m > 0:
                self._Aa_enkf, X5 = self.compute_analysis(self._Af, self._Af, z, H, R)
            else:
                X5 = None

            if self._X5 is not None: np.dot(self._X5, X5, out=self._X5)
            else: self._X5 = X5
            self._haveX5 = True

            self._Af = self._Aa_enkf # For global analysis these two coincide
            self._num_z+= m

        # Localized analysis
        else:
            if self._X5 is None:
                self._X5 = np.empty((num_domains, N, N))
                self._haveX5 = np.zeros(num_domains, dtype=np.bool)

            if m > 0:
                for d_i, d in enumerate(self._domains):
                    if d is None: continue
                    start_d = self._loc_statesize_limits[d_i][0]
                    n_d = self._loc_statesize_limits[d_i][1]

                    z_used, H_d, R_d = self.localization_strategy.get_obsdata(d, H, R)
                    m_d = len(z_used)

                    A_d = self._Aa_enkf[start_d:start_d + n_d, :]
                    assert A_d.shape == (n_d, N)

                    # Have some observations to assimilate for this local domain
                    if m_d > 0:
                        assert H_d.shape == (m_d, n)
                        assert R_d.shape == (m_d, m_d)

                        # Perform local analysis on the domain
                        z_d = z[z_used]
                        A_d_up, X5_d = self.compute_analysis(self._Af, A_d, z_d, H_d, R_d)
                        assert A_d_up.shape == (n_d, N)
                        assert X5_d.shape == (N, N)

                        if not self._haveX5[d_i]:
                            self._X5[d_i,:,:] = X5_d
                            self._haveX5[d_i] = True
                        else:
                            self._X5[d_i, :, :] = np.dot(self._X5[d_i,:,:], X5_d)

                        # Todo: Do this in-place?
                        self._Aa_enkf[start_d:start_d + n_d, :] = A_d_up

                        self._num_z[d_i] += m_d



            # After we're done with the analysis, reconstruct the global ensemble as we are going to
            # need it either at next assimmilate() call or in end_analysis().
            for d_i, d in enumerate(self._domains):
                if d is None: continue
                start_d = self._loc_statesize_limits[d_i][0]
                n_d = self._loc_statesize_limits[d_i][1]
                self.localization_strategy.put_state(d, self._Aa_enkf[start_d:start_d + n_d, :], self._Af)

        self._num_z+= m



    def end_analysis(self, on_smoother_result=None):
        n, N = self._n, self._N  # State and ensemble size
        k = len(self._smoother_data)
        num_domains = len(self._domains)

        # First the lagged smoother by updating all previous EnKF states with the
        # transformation matrix from the current analysis update. For lag=0 this
        # does nothing
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

            if num_domains == 0:
                if self._X5 is not None:
                    if self.forgetting_factor != 1.0:
                        X4 = self._X5 - eyeN
                        X4 *= self.forgetting_factor
                        self._X5 = eyeN + X4
                    Aj = Aj.dot(self._X5, out=Aj)
                if Aj_is_result: As = Aj

            # For local analyses, Aj is an SxN array with S being the sum of state vector
            # sizes for all local domains (and it may be fo that S != the global state size `n`
            # if there is some padding applied to domains!) and N is th ensemble size.
            # Likewise, X5 is an DxNxN array.
            else:
                As = np.empty((n, N)) if Aj_is_result else None

                for d_i, d in enumerate(self._domains):
                    if d is None: continue

                    start_d = self._loc_statesize_limits[d_i][0]
                    n_d = self._loc_statesize_limits[d_i][1]

                    Ajd = Aj[start_d:start_d + n_d, :]
                    assert Ajd.shape == (n_d, N)

                    if self._haveX5[d_i]: # and self._last_analysis_j[d_i] <= j:
                        X5d = self._X5[d_i, :, :]
                        if self.forgetting_factor != 1.0:
                            X4 = X5d - eyeN
                            X4 *= self.forgetting_factor
                            X5d = eyeN + X4
                            self._X5[d_i, :, :] = X5d

                        Ajd[:] = Ajd.dot(X5d)

                    if Aj_is_result: self.localization_strategy.put_state(d, Ajd, As)


            if As is not None:
                if on_smoother_result is not None:
                    on_smoother_result(ensemble.mean(As), As, tj)
                #print("Remove Aj")
                self._cache.remove(Aj_handle)


        #
        if num_domains == 0 and self._num_z > 0:
            self._last_analysis_j = k
        elif num_domains > 0:
            for d_i, d in enumerate(self._domains):
                if d is None: continue
                if self._num_z[d_i] > 0: self._last_analysis_j[d_i] = k



        # Done with smoothing, if there was any. Store analysis result for this update step
        # and return it.
        Aa_handle = self._cache.put(self._Aa_enkf)

        #print("Add Aj", k)
        self._smoother_data.append((Aa_handle, self._t))

        self._X5 = None
        self._haveX5 = None
        self._Aa_enkf = None

        # The reconstructed
        return self._Af


    def finish(self, on_smoother_result):
        n, N = self._n, self._N  # State and ensemble size
        k = len(self._smoother_data)
        num_domains = len(self._domains)

        As_buffer = np.empty((n, N)) if num_domains > 0 else None

        # First the lagged smoother by updating all previous EnKF states with the
        # transformation matrix from the current analysis update. For lag=0 this
        # does nothing
        for j in range(k - 1, k - 1 - self._lag, -1):
            if j < 0: break
            Aj_handle, tj = self._smoother_data[j]

            # print("Smoother up", k - j, Aj_handle, tj)
            Aj = self._cache.get(Aj_handle)
            assert Aj is not None

            # For global analysis both A_enkf and X5 are the global ensemble and transform
            # arrays, respectively, so we only need to compute Aj*X5
            As = None
            if num_domains == 0:
                As = Aj
            else:
                for d_i, d in enumerate(self._domains):
                    if d is None: continue
                    start_d = self._loc_statesize_limits[d_i][0]
                    n_d = self._loc_statesize_limits[d_i][1]
                    Ajd = Aj[start_d:start_d + n_d, :]
                    assert Ajd.shape == (n_d, N)
                    self.localization_strategy.put_state(d, Ajd, As_buffer)
                As = As_buffer

            if on_smoother_result is not None:
                on_smoother_result(ensemble.mean(As), As, tj)
                # print("Remove Aj")
            self._cache.remove(Aj_handle)



    def compute_analysis(self, A_global, A, z, H, R, Aa_out=None, X5_out=None):
        n, N = A.shape  # State and ensemble size
        m = H.shape[0]  # Number of observations

        #rho = 1 - (self.cov_inflation - 1)

        X5, X5s = self.variant.compute_ensemble_transform(A_global, A, z, H, R,
                                                          self.cov_inflation,
                                                          lstrategy=self.localization_strategy,
                                                          out=X5_out)
        Aa = A.dot(X5, out=Aa_out)
        return Aa, X5s if X5s is not None else X5








