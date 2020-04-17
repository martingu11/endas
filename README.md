# EnDAS

---
**The master branch is currently being reworked towards the C++/C/Python implementation.** The `0.1-py` branch
contains the original Python version. 
---

Ensemble Data Assimilation for C++, C and (eventually) Python.

## Overview

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

For installation instruction see [this page](INSTALL.md).
