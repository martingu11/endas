"""
Toy models used for testing the data assimilation algorithms.
"""

import numpy as np

__all__ = [ 'Lorenz95' ]


# Used by the Lorenz95 evolution operator to represent model trajectory.
class L95Trajectory:
    __slots__ = ['dt', 'k', 'x']


class Lorenz95:
    """
    Lorenz 95 evolution operator.

    The Lorenz 95 model is a dynamical system formulated by Edward Lorenz in [1]. It is defined as follows.
    For :math:`i = 1, ..., N`

    .. math::

        \frac {dx_{i}}{dt} = (x_{i+1}-x_{i-2})x_{i-1}-x_{i}+F

    where it is assumed that :math:`x_{-1}=x_{N-1},x_{0}=x_{N}` and :math:`x_{N+1}=x_{1}`. Here :math:`x_i`
    is the state of the system and :math:`F` is a forcing constant. :math:`F=8` is a common value known to cause
    chaotic behavior. The model is mostly useful for testing the performance of various assimilation techniques.

    References:
        [1] Lorenz, Edward (1996). "Predictability – A problem partly solved"
    """

    def __init__(self, n=40, F=8):
        self._n = n
        self._F = F

        # Precomputed array indices for the ODE expressions to avoid explicit loops later on
        self._i_m1 = np.array([i % n for i in range(-1, n-1)])  # x_i-1
        self._i_m2 = np.array([i % n for i in range(-2, n-2)])  # x_i-2
        self._i_p1 = np.array([i % n for i in range( 1, n+1)])  # x_i+1
        self._i_p2 = np.array([i % n for i in range( 2, n+2)])  # x_i+2


    def __call__(self, x, dt):
        if len(x.shape) == 1: x = x.reshape(-1, 1)
        n, N = x.shape
        assert n == self._n

        def l95(x):
            return -x[self._i_m2]*x[self._i_m1] + x[self._i_m1]*x[self._i_p1] - x + self._F

        # State propagation using the 4th-order Runge-Kutta method. Store what we need for the tl() and ad()
        # computations

        trj = L95Trajectory()
        trj.dt = dt
        trj.k = np.zeros((N, 4, n))
        trj.x = np.zeros((N, n))

        for i in range(N):
            trj.k[i,0,:] = dt * l95(x[:,i])
            trj.k[i,1,:] = dt * l95(x[:,i] + trj.k[i,0,:] / 2)
            trj.k[i,2,:] = dt * l95(x[:,i] + trj.k[i,1,:] / 2)
            trj.k[i,3,:] = dt * l95(x[:,i] + trj.k[i,2,:])
            x[:,i] += (trj.k[i,0,:] + 2*trj.k[i,1,:] + 2*trj.k[i,2,:] + trj.k[i,3,:]) / 6.0
            trj.x[i,:] = x[:,i]

        return trj


    def dot(self, trajectory, x, out=None):
        if len(x.shape) == 1: x = x.reshape(-1, 1)
        n, N = x.shape
        assert n == self._n

        assert isinstance(trajectory, L95Trajectory)

        def l95_tl(x, dx):
            return -x[self._i_m1]*dx[self._i_m2] + (x[self._i_p1]-x[self._i_m2])*dx[self._i_m1] - dx + \
                   x[self._i_m1]*dx[self._i_p1]

        if out is None:
            out = np.zeros((n, N))
        else:
            if N == 1: out = out.reshape(-1, 1)
            assert out.shape == x.shape

        for i in range(N):
            ti = i % trajectory.x.shape[0]
            dK1 = trajectory.dt * l95_tl(trajectory.x[ti,:], x[:,i])
            dK2 = trajectory.dt * l95_tl(trajectory.x[ti,:] + trajectory.k[ti,0,:] / 2, x[:,i] + dK1 / 2)
            dK3 = trajectory.dt * l95_tl(trajectory.x[ti,:] + trajectory.k[ti,1,:] / 2, x[:,i] + dK2 / 2)
            dK4 = trajectory.dt * l95_tl(trajectory.x[ti,:] + trajectory.k[ti,2,:], x[:,i] + dK3)
            out[:,i] = x[:,i] + (dK1 + 2 * dK2 + 2 * dK3 + dK4) / 6.0

        return out


    def adjdot(self, trajectory, x, out=None):
        if len(x.shape) == 1: x = x.reshape(-1, 1)
        n, N = x.shape
        assert n == self._n

        assert isinstance(trajectory, L95Trajectory)

        def l95_ad(x, dx):
            return x[self._i_m2]*dx[self._i_m1] + (x[self._i_p2]-x[self._i_m1])*dx[self._i_p1] - dx - \
                   x[self._i_p1]*dx[self._i_p2]

        def adK1(dx, i): return trajectory.dt * l95_ad(trajectory.x[i,:], dx);

        def adK2(dx, i):
            aux = l95_ad(trajectory.x[i,:] + trajectory.k[i,0,:]/2, dx)
            return trajectory.dt * (aux + adK1(aux, i) / 2)

        def adK3(dx, i):
            aux = l95_ad(trajectory.x[i,:] + trajectory.k[i,1,:]/2, dx)
            return trajectory.dt * (aux + adK2(aux, i) / 2)

        def adK4(dx, i):
            aux = l95_ad(trajectory.x[i,:] + trajectory.k[i,2,:], dx)
            return trajectory.dt * (aux + adK3(aux, i))

        if out is None:
            out = np.zeros((n, N))
        else:
            if N == 1: out = out.reshape(-1, 1)
            assert out.shape == x.shape

        for i in range(N):
            ti = i % trajectory.x.shape[0]
            dx = x[:,i]
            out[:,i] = dx + (adK1(dx, ti) + 2 * adK2(dx, ti) + 2 * adK3(dx, ti) + adK4(dx, ti)) / 6.0

        return out
