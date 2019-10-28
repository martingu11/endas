"""
Ensemble Kalman Filter.
"""

__all__ = ['EnKF']

import numpy as np
from scipy import linalg

from .enkf_base import EnKFVariant
from endas import ensemble


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




