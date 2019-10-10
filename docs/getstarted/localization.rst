.. _localization:

Localization
============

Under construction.



Truncation of covariance matrix
-------------------------------



.. _localization_domain:

Domain-based localization
-------------------------

Under this scheme, the global state is decomposed into local domains based on the (geographic)
location that each state vector element represents. Furthermore, only observations within the de-correlation length
are used for computing the analysis update of each local domain. The de-correlation length is given either as a
radius (i.e. a single scalar) or more typically via a taper function (:class:`endas.localization.TaperFn`).


