.. _da-algorithms:

DA Algorithms (endas.algorithms)
================================

Data assimilation algorithms.


.. currentmodule:: endas.algorithms

.. rubric:: Implemented algorithms

.. autosummary::
    :toctree: generated/
    :template: class.rst
    :nosignatures:

    KalmanFilter
    EnsembleKalmanFilter


.. rubric:: EnsembleKalmanFilter variants

The :class:`EnsembleKalmanFilter` class is a generic EnKF implementation that is responsible for the execution of the
sequential filter updates and for smoothing. The actual implementation of the analysis update is carried out using an
implementation of the :class:`EnKFVariant` interface class. The following EnKF variants are implemented:

.. autosummary::
    :toctree: generated/
    :template: class.rst
    :nosignatures:

    EnKF
    ESTKF





