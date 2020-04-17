"""
Extended Kalman Filter.
"""


__all__ = ['KalmanFilter']

import numpy as np
from scipy import linalg
from python.endas import CovarianceOperator
from python.endas import arraycache, cov


class KFSmootherData:
    __slots__ = ['xf', 'Pf', 'xa', 'Pa', 'trj', 't']

    def __init__(self, xf=None, Pf=None, xa=None, Pa=None, trj=None, t=None):
        self.xf = xf
        self.Pf = Pf
        self.xa = xa
        self.Pa = Pa
        self.trj = trj
        self.t = t


class KalmanFilter:
    """
    Kalman Filter and Smoother.
    """

    def __init__(self, model, model_tl, model_adj,
                 cache=None, lag=0, forgetting_factor=1.0):
        """
        Kalman Filter and Smoother.

        Args:
            model     : Callable with signature `model(x, dt) -> tladj_data` implementing forward propagation
                        of the state from time `t` to time `t+dt`. The return value `tladj_data` is
                        passed to the tangent linear and adjoint when called.
            model_tl  : Callable with signature `model(tladj_data, x)` implementing the dot product of the
                        tangent linear of the model with `x`.
            model_adj : Callable with signature `model(tladj_data, x)` implementing the dot product of the
                        adjoint of the model with `x`.
            cache:
            lag:
            forgetting_factor:
        """

        self._lag = lag
        assert lag is None or lag >= 0

        self._M = model
        self._Mtl = model_tl
        self._Madj = model_adj

        self._cache = cache if cache is not None else arraycache.ArrayCache()
        # Cache must be ArrayCache instance or subclass
        assert isinstance(self._cache, arraycache.ArrayCache)

        self.forgetting_factor = forgetting_factor
        assert forgetting_factor > 0 and forgetting_factor <= 1.0

        self._data = []
        self._nassim_this_step = 0
        self._xf_handle = None
        self._Pf_handle = None


    def smoother_begin(self, x0, P0, t0):

        if self._lag == 0: return

        P0 = cov.to_matrix(P0, force_dense=True)  # Avoid sparse matrices here for performance reasons!
        x0_handle = self._cache.put(x0)
        P0_handle = self._cache.put(P0)
        self._data.append(KFSmootherData(None, None, x0_handle, P0_handle, None, t0))


    def forecast(self, xb, Pb, Q, dt):
        """
        Implements the forecast step of the Kalman Filter.

        Args:
            xb:  Background state vector
            Pb:  Background error covariance matrix
            Q:   Model error covariance matrix. Can be `None` for perfect model
            dt:  Time increment. This is simply passed to the model

        Returns:

        """
        assert xb.ndim == 1
        n = xb.shape[0]

        Pb = cov.to_matrix(Pb, force_dense=True) # Avoid sparse matrices here for performance reasons!
        Q = cov.to_matrix(Q, force_dense=True) if Q is not None else None

        # Move state estimate and covariance from xb to xk
        trj = self._M(xb, dt)
        Pf = self._Mtl(trj, Pb)
        Pf = self._Madj(trj, Pf)
        if Q is not None: Pf += Q

        if self._lag is None or self._lag > 0:
            assert len(self._data) > 0, "forecast() called before smoother_begin()"
            self._data[-1].trj = trj

        return xb.ravel(), Pf


    def begin_analysis(self, x, P, t):
        """
        Begins the analysis update of the Kalman Filter and Smoother.
        Args:
            x: Forecast (prior) state vector
            P: Forecast (prior) error covariance matrix
            t: Analysis time

        Returns: nothing
        """

        self._xf = x
        self._Pf = cov.to_matrix(P, force_dense=True)
        self._t = t
        self._xa = self._xf
        self._Pa = self._Pf

        P = cov.to_matrix(P, force_dense=True)  # Avoid sparse matrices here for performance reasons!
        self._xf_handle = self._cache.put(x) if self._lag > 0 else None
        self._Pf_handle = self._cache.put(P) if self._lag > 0 else None


    def assimilate(self, z, H, R):
        """
        Performs assimilation of observations during the analysis update.
        Please note that `assimilate()` can only be called from within the analysis step,
        i.e. between the `begin_analysis()` and `end_analysis()` calls.

        Args:
            z:  Vector of observations.
            H:  Observation operator at time `t`.
            R:  Observation covariance matrix.

        Returns: nothing
        """

        if z is not None and len(z) > 0:
            H = H.to_matrix(force_dense=True)

            # Observation noise covariance
            F = H.dot(self._Pf).dot(H.T)

            if isinstance(R, np.ndarray): np.add(F, R, out=F)
            elif isinstance(R, CovarianceOperator): R.add_to(F)
            else:
                raise TypeError("R must be a NumPy array or endas.CovarianceOperator")


            # State update as xk + Cp*H'*F^-1*dz
            dz = z - H.dot(self._xf)
            self._xa = self._xf + self._Pf.dot(H.T).dot(linalg.solve(F, dz, sym_pos=True, overwrite_b=True))
            self._xa = self._xa.ravel()

            # Covariance estimate update as Cp - Cp*H'*F^-1*H*Cp
            self._Pa = self._Pf - self._Pf.dot(H.T).dot(linalg.solve(F, H.dot(self._Pf), sym_pos=True, overwrite_a=True, overwrite_b=True))
        else:
            self._xa = self._xf
            self._Pa = self._Pf


    def end_analysis(self, on_smoother_result=None, result_args=tuple()):
        """
        Ends the analysis update step and returns the filter analysis solution.
        Optionally, any available smoother solution is returned as well. If smoothing is enabled and
        any smoothing soluton is available at this point, `on_smoother_result` will be called (if not
        `None) as:

            on_smoother_result(xs, Ps, t)

        where `xs` and `Ps` is the state state vector and error covariance at time `t`.

        Args:
            on_smoother_result: Callable to be executed when smoothing solution is available.

        Returns: Tuple `(xa, Pa)` containing the state vector and error covariance matrix after the
        update step.
        """
        self._xa = self._xa.ravel()

        if self._lag == 0:
            if on_smoother_result is not None: on_smoother_result(self._xa, self._Pa, self._t, result_args)
        else:
            xa_handle = self._cache.put(self._xa)
            Pa_handle = self._cache.put(self._Pa)
            self._data.append(KFSmootherData(self._xf_handle, self._Pf_handle, xa_handle, Pa_handle, None, self._t))

        return self._xa, self._Pa


    def smoother_finish(self, on_smoother_result, result_args=tuple()):
        """
        Finalizes any pending smoothing solutions and calls `on_smoother_result`.
        For each available smoother solution, `on_smoother_result` will be called as:

            on_smoother_result(xs, Ps, t)

        where `xs` and `Ps` is the state state vector and error covariance at time `t`.

        Args:
            on_smoother_result: Callable to be executed for all available smoothing solutions.

        Returns: nothing
        """

        if self._lag == 0: return

        smd = self._data[-1]
        #_, _, xa_handle, Pa_handle, trj, t = self._data[-1]

        # Todo: implement lock-release into ArrayCache!
        xs = self._cache.get(smd.xa)
        Ps = self._cache.get(smd.Pa)
        self._cache.remove(smd.xa)
        self._cache.remove(smd.Pa)

        on_smoother_result(xs, Ps, smd.t, result_args)

        for k in range(len(self._data)-2, -1, -1):
            smdk1 = self._data[k+1]
            smdk = self._data[k]
            #xf_handle, Pf_handle = self._data[k+1][0:2]
            #xa_handle, Pa_handle, trj, t = self._data[k][2:]

            xf = self._cache.get(smdk1.xf)
            Pf = self._cache.get(smdk1.Pf)
            xa = self._cache.get(smdk.xa)
            Pa = self._cache.get(smdk.Pa)

            J = linalg.solve(a=Pf, b=self._Mtl(smdk.trj, Pa), assume_a='sym').T

            J*= self.forgetting_factor

            xs = xa + J.dot(xs - xf)
            Ps = Pa + J.dot(J.dot(Ps - Pf).T)

            self._cache.remove(smdk1.xf)
            self._cache.remove(smdk1.Pf)
            self._cache.remove(smdk.xa)
            self._cache.remove(smdk.Pa)

            on_smoother_result(xs, Ps, smdk.t, result_args)


