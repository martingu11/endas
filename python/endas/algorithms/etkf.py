"""
Ensemble Transform Kalman Filter (ETKF)
=======================================


"""

__all__ = ['ESTKF']

import math
import numpy as np
from scipy import linalg

from .enkf import EnKFVariant
from .. import ensemble


class ESTKF(EnKFVariant):
    """
    Error Subspace Transform Kalman Filter.

    """

    def __init__(self, rotation=False):
        self.rotation=rotation


    def process_global_ensemble(self, Ag, H):
        xg = ensemble.mean(Ag)
        Hx = H.dot(xg)
        HA = H.dot(Ag)
        return Hx, HA


    def ensemble_transform(self, A, z, H, R, Ag_data, inflation, lstrategy, out=None):
        n, N = A.shape  # State and ensemble size
        #m = H.shape[0]  # Number of observations

        rho = 1.0 - (inflation - 1.0)
        #rho_s = 1.0 - (rho - 1.0)

        Hx, HA = Ag_data

        a = (1.0 / N) * (1.0 / (1.0 / math.sqrt(N) + 1))
        T = np.full((N, N - 1), -a)
        np.fill_diagonal(T, 1.0 - a)
        T[-1, :] = -1.0 / math.sqrt(N)

        HL = HA.dot(T)
        RinvHL = R.solve(HL)

        Ainv = HL.T.dot(RinvHL)
        np.fill_diagonal(Ainv, Ainv.diagonal() + (rho * (N - 1)))
        # np.fill_diagonal(Ainv, Ainv.diagonal() + ((m - 1)))

        dz = z - Hx
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
