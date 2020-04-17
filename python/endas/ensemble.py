"""
Basic utility functions for working with ensembles.
"""

import math
import numpy as np
from scipy import linalg, sparse


def mean(A):
    # Returns ensemble mean of the ensemble `A`
    return np.mean(A, axis=1)


def to_anomaly(E, Eu=None, normalize=False, out=None):
    """
    Computes ensemble anomaly from the ensemble and its mean. Optionally the anomalies can
    be normalized by ``sqrt(m-1)`` where `m` is the ensemble size.

    Args:
        E : NxM array of M ensemble members (i.e. each column corresponds to ensemble member
            state vector of length N).
        Eu : Array of length M representing the ensemble mean
        normalize : If True, the anomalies are divided by ``sqrt(m-1)`` where `m` is the ensemble size
        out : Array of shape NxM where the result is stored or None if new array should be allocated

    Returns: NxM array of anomalies
    """
    n, m = E.shape
    if Eu is None: Eu = np.mean(E, axis=1)

    A = np.subtract(E, Eu.reshape(n, 1), out=out)
    if normalize: np.divide(A, math.sqrt(m - 1), out=A)
    return A



def from_anomaly(E, Eu, out=None):
    """
    Computes ensemble from ensemble mean and anomaly vectors.
    This is an inverse of to_anomaly().
    """
    #Todo: Add `normalize`
    n, m = E.shape
    return np.add(E, Eu.reshape(n, 1), out=out)



def center(E, out=None):
    """
    Adjusts the ensemble to zero mean.
    """
    u = np.mean(E, axis=1).reshape(E.shape[0], 1)
    return np.subtract(E, u, out=out)


def standardize(E, out=None):
    """
    Adjusts the ensemble to zero mean and standard deviation of 1.
    """
    u = np.mean(E, axis=1).reshape(-1, 1)
    E = np.subtract(E, u, out=out)

    sd = np.std(E, axis=1).reshape(-1, 1)
    E = np.divide(E, sd, out=out)
    return E


def inflate(A, inflation,  Au = None, out=None):

    if Au is None: Au = mean(A)
    AX = to_anomaly(A, Eu=Au, out=out)
    AX = np.multiply(AX, inflation, out=AX)
    A = from_anomaly(AX, Eu=Au, out=AX)
    return A



def generate(N, u, cov):
    """
    Generates new ensemble with given covariance.

    Args:
        N : Size of the ensemble to generate
        u : Desired ensemble mean as a scalar or vector
        cov : Covariance of the ensemble to generate as CovarianceMatrix instance
        centered : Ensure the ensemble has zero mean

    Returns:
        New nxN array where `n` is the state size (``n = cov.shape[0]``) and `N` is the
        requested ensemble size
    """
    E = cov.random_multivariate_normal(N)
    if isinstance(u, np.ndarray): u = u.reshape(-1, 1)
    np.add(E, u, out=E)
    return E




def optimal_resample(E, alpha=1):
    """
    Resamples an augmented ensemble of size n x alpha*N to the desired size n x N
    so that the rank of the resulting ensemble is (or is close to)
    rank(E') = min(n, N).

    To use this, generate a larger ensemble, (i.e. for alpha=6 generate an ensemble
    of size 6N members, therefore an n x 6N matrix). Then, call optimal_resample()
    to obtain the NxN ensemble with as high rank as possible.
    """

    n, Nalpha = E.shape
    N = int(Nalpha // alpha)
    assert N*alpha == Nalpha # Input ensemble size must be a multiple of alpha!
    if N == 1 or n < Nalpha: return E

    if n > Nalpha:
        U, sdiag, Vt = linalg.svd(E, full_matrices=False, overwrite_a=True)
        U = U[:, 0:N]
        sdiag = sdiag[0:N]
        np.multiply(sdiag, math.sqrt(1.0 / alpha), out=sdiag)

        Y = np.random.standard_normal(size=(N, N))
        Q, R = linalg.qr(Y, overwrite_a=True)

        Rdiag = np.diagonal(R, 0)
        Xdiag = Rdiag / np.fabs(Rdiag)
        X = np.diag(np.reciprocal(Xdiag), 0)
        Q = Q.dot(X)
        Q = Q[0:N, :]

        # Note: This implements S Q^T where S is diagonal matrix made from `sdiag`.
        # Avoiding use of scipy.sparse here as it is slow for this!
        #sQT = np.multiply(Q.T, sdiag.reshape(-1, 1))

        s = sparse.diags(sdiag)
        #E2 = U.dot(sQT)
        E2 = U.dot(s.dot(Q.T))
        return E2

    # Todo: Should do something about this case also!
    else:
        return E[:,0:N]



