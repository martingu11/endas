Localization (endas.localization)
=================================

Localization of the analysis update.

This module implements various forms of localization for ensemble-based data assimilation methods.

.. currentmodule:: endas.localization

Classes and functions
----------------------------------------

.. rubric:: Interface classes for localization

.. autosummary::
    :toctree: generated/
    :nosignatures:
    :template: class.rst

    CoordinateSystem
    TaperFn
    SpatialQuery
    StateSpacePartitioning
    DomainLocalization


.. rubric:: Coordinate systems

.. autosummary::
    :toctree: generated/
    :template: class.rst
    :nosignatures:

    cs.EuclideanCS
    cs.LatLonCS


.. rubric:: Covariance tapering functions

.. autosummary::
    :toctree: generated/
    :template: class.rst
    :nosignatures:

    taper.GaspariCohn
    taper.Linear
    taper.Spherical


.. rubric:: Classes implementing domain localization


Instead, this is delegated to classes that implement the :class:`StateSpacePartitioning` interface. Currently, EnDAS
contains the following state space partitioning schemes:

.. autosummary::
    :toctree: generated/
    :template: class.rst
    :nosignatures:

    GenericStateSpace1d
    grid.Grid2d


Implementing custom state space geometry
----------------------------------------


Implementing custom localization scheme
---------------------------------------
