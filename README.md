# EnDAS

Ensemble Data Assimilation for Python 3.x.

## Overview

EnDAS is a data assimilation library for Python that focuses on ensemble data assimilation algorithms 
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

## Requirements

EnDAS is written in Python 3, using NumPy and SciPy. Therefore, if you cannot use either, EnDAS is unfortunately
not for you. Some parts of EnDAS are implemented as C/C++ extensions (using Cython) so you will need working 
C/C++ compiler and, if building from the GitHub source, Cython installed. EnDAS will work on any Python interpreter 
implementation compatible with the Cython such as CPython or PyPy. For more information see the 
[installation instructions](INSTALL.md).  


## Documentation

Documentation can be found at https://endas.readthedocs.io/en/latest.

Code examples can be found in the `examples` directory of the EnDAS repository on GitHub:
[https://github.com/martingu11/endas]([https://github.com/martingu11/endas]).

## Installation

For installation instruction see [this page](INSTALL.md).
