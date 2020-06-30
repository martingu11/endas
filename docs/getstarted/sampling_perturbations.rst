.. _sampling-perturbations:

Sampling of initial ensemble and perturbations
==============================================


Using covariance matrix on small problems
-----------------------------------------

When sampling from a known covariance matrix :math:`\boldsymbol{\Sigma}`, the :class:`endas::DenseCovarianceOp` 
class can be conveniently used to draw realizations from :math:`\mathcal{N}(\mathbf{0}, \boldsymbol{\Sigma})` 
using the ``randomMultivariateNormal()`` method. 

.. rubric:: C++

.. code-block:: cpp
   
   // The covariance matrix
   endas::Matrix cov = ...; 

   endas::DenseCovarianceOp covOp(cov);

   // Sample ensemble perturbations into an array of size nxN (n=state size, N=ensemble size)
   Array2d E(n, N);  
   covOp.randomMultivariateNormal(E);

   // Add desired mean
   E += 1.2345;


.. note:: 
   This approach is only feasible for small problems due to the size of covariance matrix (:math:`n \times n`,
   where *n* is the size of the state space) and computational costs of the Cholesky decomposition and matrix 
   products involved.


The :class:`~endas::DenseCovarianceOp` class can also be initialized from a covariance function 
(:class:`~endas::CovarianceFn`). In this case, the covariance matrix is created internally using the covariance 
information provided by the function.



Generating random Gaussian fields
---------------------------------

For larger state spaces, covariance matrices cannot be used directly. If the state space is organized on a 
two-dimensional grid, random Gaussian fields can be generated according to a covariance function via the
:class:`~endas::GaussianRandomField` class. GaussianRandomField currently only supports isotropic covariance
functions such that the covariance only depends on distance, not direction.

The "spatial" EnDAS module contains implementation of several widely used covariance functions such as 
:any:`~endas::ExponentialCovFn`, :any:`~endas::GaussianCovFn` and :any:`~endas::SphericalCovFn`.


.. rubric:: C++

.. code-block:: cpp
   
   // Gaussian covariance function will correlation length 10 and variance 1
   auto covFn = endas::GaussianCovFn(10, 1); 

   // Gaussian random field sampler over a 100x100 grid
   endas::GaussianRandomField grf(100, 100, covFn);

   // Sample one realization
   Array2d E(100, 100);  
   grf.sample(E);

   // Add desired mean
   E += 1.2345;


:class:`~endas::GaussianRandomField` uses Fast Fourier Transform to generate the perturbations and
two independent realizations are obtained simultaneously (one in the real, one in the imaginary part). 
It is therefore possible to sample two realizations at once at the same cost as sampling a single one
by supplying two arrays to ``sample()``. For convenience, the :func:`~endas::GaussianRandomField::sampleEnsemble()` 
method can be used to draw many realizations in a single call.


Exact second-order sampling
---------------------------

See :class:`~endas::SecondOrderExactSample`.
