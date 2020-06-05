Simple example
==============

Example below uses a simple 3-dimensional system to demonstrate data assimilation with
EnDAS. The system evolves according to a linear model

.. math::

   \mathbf{M} = 
       \begin{bmatrix}
           1  &  0                                 & 0 \\ 
           0  &   \cos \left(\frac{\pi}{6} \right) & \sin \left(\frac{\pi}{6} \right) \\
           0  &  -\sin \left(\frac{\pi}{6} \right) & \cos \left(\frac{\pi}{6} \right)
       \end{bmatrix}

which is a simple rotation around the first state dimension by :math:`\pi / 6`. Therefore, 
the system is periodic with a period of 12 steps. The data assimilation start from the state 
:math:`(x_1, x_2, x_3)^T = (0, 0, 1)^{\mathrm{T}}` and observations are assimilated on every third time 
step. The observation operator is given by matrix 

.. math::

   \mathbf{H} = \begin{bmatrix}1  &  1 & 0 \end{bmatrix}

so that the quantity :math:`x_1 + x_2` is observed while the state variable :math:`x_3` remains
unobserved.


.. figure:: /images/example_simple.png
   
   Data fusion of the simple rotating system. Black line represents the "true" system state (:math:`x_1+x_2`), 
   purple dots are the assimilated observations. The state estimate and uncertainty (standard deviation) are 
   shown in green.


Below is the full source code:

.. literalinclude:: /../examples/SimpleExample.cpp
   :language: cpp

