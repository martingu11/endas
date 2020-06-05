EnDAS
=====

Ensemble Data Assimilation for C++, C and (eventually) Python.

Overview
--------

EnDAS is a data assimilation library that focuses on ensemble data assimilation algorithms 
(although few others are also included for comparison). The main features are:

-  Ensemble Kalman Filters and Smoothers including
   - Traditional/stochastic EnKF
   - Square root EnKF
   - Ensemble Transform Kalman Filters
   - Variational Ensemble Kalman Filter (still experimental)  
-  Traditional Kalman Filter and Smoother
-  Distance-based localization of the analysis update, including few popular covariance tapering functions
-  Non-intrusive filtering (and smoothing) API
-  Utilities for generating random fields





Requirements
------------

EnDAS is written in C++11 and its only required dependency is `Eigen <http://eigen.tuxfamily.org>`_. 
Please see the :doc:`Installation instructions page <installation>` for more information.


Installation
------------

See the :doc:`Installation instructions page <installation>`.


Where to start?
---------------



Code examples can be found in the **examples** directory of the
`EnDAS repository on GitHub <https://github.com/martingu11/endas>`_.




Indices and tables
==================

* :ref:`genindex`
* :ref:`search`



.. toctree::
   :maxdepth: 1
   :hidden:

   self
   installation
   getstarted/index
   examples/index
   cppapi/cppapi_root
