
EnDAS
=====

Ensemble Data Assimilation for Python 3.x.

Overview
--------

EnDAS is a data assimilation library for Python that focuses on ensemble data assimilation algorithms
(although few others are also included for comparison). The main features are:

*  Ensemble Kalman Filters and Smoothers including

     * Traditional/stochastic EnKF
     * Square root EnKF
     * Ensemble Transform Kalman Filters
     * Variational Ensemble Kalman Filter (still experimental)

*  Traditional Kalman Filter and Smoother
*  Distance--based localization of the analysis update, including few popular covariance tapering functions
*  Non-intrusive filtering (and smoothing) API
*  Utilities for generating random fields

Requirements
------------

EnDAS is written in Python 3, using NumPy and SciPy. Therefore, if you cannot use either, EnDAS is unfortunately
not for you. Some parts of EnDAS are implemented in Cython, these are however optional and pure Python/NumPy
implementations also exist as a fall-back. For more information see the
:doc:`Installation instructions page <installation>`.

Installation
------------

For installation instruction see the :doc:`Installation instructions page <installation>`.


Documentation
-------------

TO BE COMPLETETED.

Code examples can be found in the **examples** directory of the
`EnDAS repository on GitHub <https://github.com/martingu11/endas>`_.





.. toctree::
   :maxdepth: 1
   :hidden:

   self
   installation
   api_reference





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
