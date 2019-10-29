"""
Example models
==============

Forward propagation models used in the EnDAS examples.

Please note that none of the models are tuned much for performance, they are included for
demonstration purposes only.
"""

from abc import ABCMeta, abstractmethod

__all__ = ['ForwardModel']


class ForwardModel(metaclass=ABCMeta):
    """
    Abstract base class for the test models included with EnDAS.

    Please note that you do not need to derive from this abstract base when implementing your
    own model, it only exists to ease the implementation of the example models. In fact, the
    forward propagation of the state does not require a "model class" at all (unless you want
    to use the Kalman Filter algorithm). EnDAS does not care how the ensemble is propagated
    forward.

    The model implements forward propagation via the function call operator and optionally
    the tangent linear and adjoint via the `dot()` and `adjdot()` methods.
    """

    @property
    @abstractmethod
    def ndim(self):
        """
        The number of dimensions of the model state
        """
        pass

    @abstractmethod
    def __call__(self, x, dt):
        """
        Propagates model from the current state to a new state at `t+dt`.

        The current model state `x` is assumed to be at the time `t` and is
        propagated to the time `t+dt`. The state is to be updated in-place, replacing
        the contents of `x` with the updated state vector.
        For linear operators or nonlinear operators that do not implement linearization, the
        method will usually return ``None``. For linearized operators (see
        :class:`LinearizedEvolutionOperator`), the method must return trajectory data for
        the tangent linear and adjoint operators. The trajectory data is specific to the
        operator implementation.

        Args:
            x  (NxM array) : One or more state vectors stored in columns. Outside of ensemble
                             filtering methods, M is likely to be 1 (only a single state vector
                             to be propagated).
            dt (float)     : Propagation time step

        Returns:
            Model trajectory at `t+dt` or ``None``.


        Subclasses must implement this method at least for M = 1.
        The default implementation handles the case of M > 1 by calling this method for each
        column of `x` individually and can therefore be used form the subclass as a convenience.
        The individual (sub)model trajectories are then joined together and returned as a
        CompositeTrajectory instance.
        """
        n, N = x.shape

        # Fallback implementation of ensemble state propagation via simple loop.
        if N > 1:
            trajectory = [ ]
            for k in range(N):
                trj = self(x[:, k], dt)
                trajectory.append(trj)

            return trajectory
        else:
            raise NotImplementedError()


    @abstractmethod
    def dot(self, trajectory, x):
        """
        Implements tangent-linear operator for the evolution model.

        Args:
          trajectory    : Trajectory at which the linearization occurs. This is the instance
                          returned from :func:`da.EvolutionOperator.__call__()`.
          x (NxM array) : Vector on which the tangent-linear operates.

        Returns:
          NxM array     : The product :math:`\mathbf{A}x`, where :math:`\mathbf{A}` is the
                          tangent-linear operator.


        Subclasses must implement this method at least for M=1 (M is the number of columns in `x`).
        The default implementation handles the case of M > 1 by calling this method for each
        column of `x` individually and can therefore be used form the subclass as a convenience.
        """

        n, N = x.shape

        if N > 1:
            assert isinstance(trajectory, list)
            assert N == len(trajectory)

            for k in range(N):
                x[:, k] = self.dot(trajectory[k], x[:, k])

            return x
        else:
            raise NotImplementedError()


    @abstractmethod
    def adjdot(self, trajectory, x):
        """
        Implements adjoint operator for the evolution model.

        Args:
            trajectory    : Trajectory at which the linearization occurs. This is the instance
                            returned from :func:`da.EvolutionOperator.__call__()`.
            x (NxM array) : Vector on which the adjoint operates.

        Returns:
            NxM array     : The product :math:`\mathbf{A}^\mathsf{T} x` where :math:`\mathbf{A}^\mathsf{T}`
                            is the adjoint operator.

        Subclasses must implement this method at least for M=1 (M is the number of columns in `x`).
        The default implementation handles the case of M > 1 by calling this method for each
        column of `x` individually and can therefore be used form the subclass as a convenience.
        """

        n, N = x.shape

        if N > 1:
            assert isinstance(trajectory, list)
            assert N == len(trajectory)

            for k in range(N):
                x[:, k] = self.adjdot(trajectory[k], x[:, k])

            return x
        else:
            raise NotImplementedError()
