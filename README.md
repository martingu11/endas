# EnDAS

Ensemble Data Assimilation for C++, C and (eventually) Python.

> **NOTE**: The master branch contains the latest C++ implementation of the library. The original 
  Python-only version is on the `0.1-py` branch.


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


## Installation

See [this page](https://endas.readthedocs.io/en/latest/installation.html) for installation instructions.

## Documentation

EnDAS documentation can be found [here](https://endas.readthedocs.io/en/latest/).

