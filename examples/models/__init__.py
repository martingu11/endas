"""
Example models
==============

Please note that none of the models are tuned much for performance, they are included for demonstration purposes only.
"""

from abc import ABC, abstractmethod

import numpy as np



class CompositeTrajectory(list):
  """
  Class used by the EvolutionOperator implementation to store multiple individual
  model trajectories. Derives from list and is meant to just be a collection of
  trajectory instances.
  """
  pass



class EvolutionOperator(ABC):
  """
  Abstract base class for generic evolution operators.

  This class represents evolution operators that only implement forward
  propagation in time and do not support any form of linearization.
  """

  @property
  @abstractmethod
  def ndim(self):
    """
    The number of dimensions of the model state
    """
    pass


  @abstractmethod
  def __call__(self, x, dt, trajectoryout=None):
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
      dt (float)     : Propagation time step.

    Returns:
      Model trajectory at `t+dt` or ``None``.


    Subclasses must implement this method at least for M = 1.
    The default implementation handles the case of M > 1 by calling this method for each
    column of `x` individually and can therefore be used form the subclass as a convenience.
    The individual (sub)model trajectories are then joined together and returned as a
    CompositeTrajectory instance.
    """

    N, M = x.shape

    # Fallback implementation of ensemble state propagation via simple loop.
    # Subclasses may (should) do better than this
    if M > 1:

      T = CompositeTrajectory() if isinstance(self, LinearizedEvolutionOperator) else None
      for m in range(M):
        trj = self(x[:,m], dt)
        if T: T.append(trj)

      return T
    else:
      raise Exception("EvolutionOperator.__call__() not implemented for M=1")




class LinearizedEvolutionOperator(EvolutionOperator):
  """
  Abstract base class for evolution operators that implement both the tangent linear
  and adjoint. The class is an extension of the EvolutionOperator abstract base.
  """

  @abstractmethod
  def dot(self, trajectory, x, out=None):
    """
    Implements tangent-linear operator for the evolution model.

    Args:
      trajectory    : Trajectory at which the linearization occurs. This is the instance
                      returned from :func:`da.EvolutionOperator.__call__()`.
      x (NxM array) : Vector on which the tangent-linear operates.

    Returns:
      NxM array     : The product :math:`\mathbf{A}x`, where :math:`\mathbf{A}` is the tangent-linear operator.


    Subclasses must implement this method at least for M=1 (M is the number of columns in `x`).
    The default implementation handles the case of M > 1 by calling this method for each
    column of `x` individually and can therefore be used form the subclass as a convenience.

    """

    N, M = x.shape

    if M > 1:
      haveCompositeTrj = isinstance(trajectory, CompositeTrajectory)
      if haveCompositeTrj: assert M == len(trajectory)

      for m in range(M):
        x[:,m] = self.dot(trajectory[m] if haveCompositeTrj else trajectory, x[:,m])

      return x
    else:
      raise Exception("LinearizedEvolutionOperator.tl() not implemented")



  @abstractmethod
  def adjdot(self, trajectory, x, out=None):
    """
    Implements adjoint operator for the evolution model.

    Args:
      trajectory    : Trajectory at which the linearization occurs. This is the instance
                      returned from :func:`da.EvolutionOperator.__call__()`.
      x (NxM array) : Vector on which the adjoint operates.

    Returns:
      NxM array     : The product :math:`\mathbf{A}^\mathsf{T} x` where :math:`\mathbf{A}^\mathsf{T}` is the adjoint operator.

    Subclasses must implement this method at least for M=1 (M is the number of columns in `x`).
    The default implementation handles the case of M > 1 by calling this method for each
    column of `x` individually and can therefore be used form the subclass as a convenience.
    """

    N, M = x.shape

    if M > 1:
      haveCompositeTrj = isinstance(trajectory, CompositeTrajectory)
      if haveCompositeTrj: assert M == len(trajectory)

      for m in range(M):
        x[:,m] = self.adjdot(trajectory[m] if haveCompositeTrj else trajectory, x[:,m])

      return x
    else:
      raise Exception("LinearizedEvolutionOperator.ad() not implemented")



