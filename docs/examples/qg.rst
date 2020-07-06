Example with the Quasi-geostrophic model
========================================

Example below demonstrates data assimilation using a 1.5 -layer Quasi-Geostrophic model.
For more information about the model and the data assimilation set up, see 
 
SAKOV, P. and OKE, P.R. (2008), A deterministic formulation of the ensemble Kalman filter: 
an alternative to ensemble square root filters. Tellus A, 60: 361-371. 
doi:10.1111/j.1600-0870.2007.00299.x.

.. figure:: /images/example_qg.jpg
   
   Data fusion of the simple rotating system. Black line represents the "true" system state (:math:`x_1+x_2`), 
   purple dots are the assimilated observations. The state estimate and uncertainty (standard deviation) are 
   shown in green.


Below is the full source code:

.. literalinclude:: /../examples/QGExample.cpp
   :language: cpp

