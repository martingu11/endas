"""
Ensemble Kalman Filter.
"""

__all__ = [ 'VariantEnKF' ]

import numpy as np
from scipy import linalg

from .enkf_base import EnKFVariant
from endas import ensemble


class VariantEnKF(EnKFVariant):
    """
    Classic (stochastic) Ensemble Kalman Filter with perturbed observations.
    """

    def compute_ensemble_transform(self, A_global, A, z, H, R, inflation, lstrategy, out=None):
        n, N = A.shape  # State and ensemble size
        m = H.shape[0]  # Number of observations

        AgX = ensemble.to_anomaly(A_global)
        HAX = H.dot(AgX)


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
        D = np.subtract(D, H.dot(A_global), out=D)
        X4 = K.dot(D, out=out)

        X5 = np.eye(N) + X4
        #np.fill_diagonal(X4, X4.diagonal() + 1.0)  # X4 + I in-place
        return X5, None




