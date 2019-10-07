"""
Toy models used for testing the data assimilation algorithms.
"""

import numpy as np
from . import ForwardModel


__all__ = [ 'Lorenz95' ]


# Used by the Lorenz95 evolution operator to represent model trajectory.
class L95Trajectory:
    __slots__ = ['dt', 'k1', 'k2', 'k3', 'k4', 'x']


class Lorenz95(ForwardModel):
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

    def __init__(self, N=40, F=8):
        self._N = N
        self._F = F

        # Precomputed array indices for the ODE expressions to avoid explicit
        # loops later on
        self._i_m1 = np.array([i % N for i in range(-1, N-1)])  # x_i-1
        self._i_m2 = np.array([i % N for i in range(-2, N-2)])  # x_i-2
        self._i_p1 = np.array([i % N for i in range( 1, N+1)])  # x_i+1
        self._i_p2 = np.array([i % N for i in range( 2, N+2)])  # x_i+2

    @property
    def ndim(self): return self._N


    def __call__(self, x, dt):
        N, M = x.shape if len(x.shape) > 1 else (len(x), 1)
        assert N == self.ndim

        # Let the base class handle the loop over ensemble members
        if M > 1: return super().__call__(x, dt)

        def l95(x):
            return -x[self._i_m2]*x[self._i_m1] + x[self._i_m1]*x[self._i_p1] - x + self._F

        # State propagation using the 4th-order Runge-Kutta method. Store what we need
        # for the tl() and ad() computations

        trj = L95Trajectory()

        trj.dt = dt
        trj.k1 = dt * l95(x)
        trj.k2 = dt * l95(x + trj.k1/2)
        trj.k3 = dt * l95(x + trj.k2/2)
        trj.k4 = dt * l95(x + trj.k3)
        x+= (trj.k1 + 2*trj.k2 + 2*trj.k3 + trj.k4) / 6.0
        trj.x = np.copy(x)

        return trj


    def dot(self, trajectory, x, out=None):

        N, M = x.shape if len(x.shape) > 1 else (len(x), 1)
        assert N == self.ndim

        # Let the base class handle the loop over ensemble members
        if M > 1: return super().dot(trajectory, x)
        assert isinstance(trajectory, L95Trajectory)


        def l95_tl(x, dx):
            return -x[self._i_m1]*dx[self._i_m2] + (x[self._i_p1]-x[self._i_m2])*dx[self._i_m1] - dx + \
                   x[self._i_m1]*dx[self._i_p1]

        dK1 = trajectory.dt * l95_tl(trajectory.x, x)
        dK2 = trajectory.dt * l95_tl(trajectory.x + trajectory.k1/2, x + dK1/2)
        dK3 = trajectory.dt * l95_tl(trajectory.x + trajectory.k2/2, x + dK2/2)
        dK4 = trajectory.dt * l95_tl(trajectory.x + trajectory.k3, x + dK3)

        return x + (dK1 + 2 * dK2 + 2 * dK3 + dK4) / 6.0



    def adjdot(self, trajectory, x, out=None):
        N, M = x.shape if len(x.shape) > 1 else (len(x), 1)
        assert N == self.ndim

        # Let the base class handle the loop over ensemble members
        if M > 1: return super().adjdot(trajectory, x, out=out)
        assert isinstance(trajectory, L95Trajectory)

        def l95_ad(x, dx):
            return x[self._i_m2]*dx[self._i_m1] + (x[self._i_p2]-x[self._i_m1])*dx[self._i_p1] - dx - \
                   x[self._i_p1]*dx[self._i_p2]

        def adK1(dx): return trajectory.dt * l95_ad(trajectory.x, dx);

        def adK2(dx):
            aux = l95_ad(trajectory.x + trajectory.k1/2, dx)
            return trajectory.dt * (aux + adK1(aux)/2)

        def adK3(dx):
            aux = l95_ad(trajectory.x  + trajectory.k2/2, dx)
            return trajectory.dt * (aux + adK2(aux)/2)

        def adK4(dx):
            aux = l95_ad(trajectory.x + trajectory.k3, dx)
            return trajectory.dt * (aux + adK3(aux))

        dx2 = x + (adK1(x) + 2 * adK2(x) + 2 * adK3(x) + adK4(x)) / 6.0
        return dx2
