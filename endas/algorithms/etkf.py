"""
Ensemble Transform Kalman Filter (ETKF)
=======================================


"""

__all__ = [ 'VariantESTKF' ]

import math
import numpy as np
from scipy import linalg

from .enkf_base import EnKFVariant
from endas import ensemble


class VariantESTKF(EnKFVariant):
    """
    Error Subspace Transform Kalman Filter.


    """

    def __init__(self, rotation=True):
        self.rotation=rotation


    def compute_ensemble_transform(self, A_global, A, z, H, R, inflation, lstrategy, out=None):
        n, N = A.shape  # State and ensemble size
        #m = H.shape[0]  # Number of observations

        rho = 1.0 - (inflation - 1.0)
        #rho_s = 1.0 - (rho - 1.0)

        x_global = ensemble.mean(A_global)

        #AgX = ensemble.to_anomaly(A_global)
        #HAX = H.dot(AgX)

        a = (1.0 / N) * (1.0 / (1.0 / math.sqrt(N) + 1))
        T = np.full((N, N - 1), -a)
        np.fill_diagonal(T, 1.0 - a)
        T[-1, :] = -1.0 / math.sqrt(N)

        HL = H.dot(A_global).dot(T)

        RinvHL = R.solve(HL)

        Ainv = HL.T.dot(RinvHL)
        np.fill_diagonal(Ainv, Ainv.diagonal() + (rho * (N - 1)))
        # np.fill_diagonal(Ainv, Ainv.diagonal() + ((m - 1)))

        dz = z - H.dot(x_global)
        w = HL.T.dot(R.solve(dz))
        w = linalg.solve(Ainv, w, overwrite_b=True)

        U, s, _ = linalg.svd(Ainv, overwrite_a=True)
        np.power(s, -0.5, out=s)
        C = U.dot(np.diag(s).dot(U.T))

        W = math.sqrt(N-1) * (C.dot(T.T))

        if self.rotation:
            Y = np.random.standard_normal(size=(N, N))
            Q, RR = linalg.qr(Y, overwrite_a=True)

            RRdiag = np.diagonal(RR, 0)
            Ediag = RRdiag / np.fabs(RRdiag)
            E = np.diag(np.reciprocal(Ediag), 0)
            Q = Q.dot(E)

            W = W.dot(Q)

        dW = w.reshape(-1, 1) + W
        G = T.dot(dW)
        G += (1.0 / N)

        Gs = (rho*T).dot(dW)
        Gs += (1.0 / N)

        return G, Gs






