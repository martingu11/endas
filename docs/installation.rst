Installation
============

EnDAS currently does not provide precompiled binaries for any platform. It is however easy to build 
from the source code, which is available from the `EnDAS GitHub repository <https://github.com/martingu11/endas>`_.


Requirements
------------

EnDAS is written in C++11 and its only required dependency is Eigen `Eigen <http://eigen.tuxfamily.org>`_. 
EnDAS can however utilize other libraries such as BLAS and LAPACK to improve performance. Their use is 
**highly reccommended** when working with large problems. CMake is used to automatically find dependencies, 
configure the source tree and generate native projects for the development environment/compiler of your choice. 

Here are the requirements in a nutshell:

**REQUIRED**:

- `CMake <https://cmake.org/>`_ version 3.13 or newer 
- `Eigen <http://eigen.tuxfamily.org>`_ version 3.3 or newer 
- Working C++ compiler with full C++11 compliance (installed on Unix-like systems by default)

Optional:

- BLAS and LAPACK (see below)
- OpenMP 
- Python 2.x development SDK (for plotting functionality in examples and benchmarks)


Eigen
^^^^^

Eigen is a header-only library that can be downloaded from the `Eigen website <http://eigen.tuxfamily.org>` or
the `GitLab repo <https://gitlab.com/libeigen/eigen>`_.


BLAS and LAPACK
^^^^^^^^^^^^^^^

BLAS and LAPACK provide highly optimized and parallelized routines for linear algebra, on which EnDAS relies 
heavily. Without BLAS/LAPACK, EnDAS will fall back on the C++ implementation provided by Eigen. The C++ 
implementation may be entirely sufficient for smaller problems and it should be noted that BLAS/LAPACK will
in any case only be used if matrices are large enough. For larger problems, the use of BLAS and LAPACK is 
reccommended.

.. note::
   LAPACK is written in Fortran and Eigen will also need the C interface called LAPACKE. Therefore, LAPACKE 
   must also be installed for EnDAS to use LAPACK.

On Linux systems that come with ``apt-get``, the `OpenBLAS implementation <https://www.openblas.net/>`_ of 
BLAS can usually be installed by::

    > sudo apt-get install libopenblas-dev

OpenBLAS provides precompiled binaries for Windows.

Similarly, the `Netlib implementation <http://www.netlib.org/lapack/index.html>`_ of LAPACK and LAPACKE can 
usually be installed by::

    > sudo apt-get install liblapack-dev liblapacke-dev


Intel MKL
^^^^^^^^^

To be written...



Building and installing from source on Linux
--------------------------------------------

The source code is hosted on GitHub: `<https://github.com/martingu11/endas>`_. After cloning the repository, you can
build and install EnDAS as follows:

Getting ready
^^^^^^^^^^^^^

If you have not yet done so, configure, build and install the Eigen library. To do so, change (``cd``) to 
the Eigen source directory, and create a ``build`` subdirectory::

    cd <eigen-dir>
    mkdir build

Change (``cd``) to the build directory and run CMake. You may want to select where the header files will be 
installed via the ``CMAKE_INSTALL_PREFIX`` option::

    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<path-to-some-folder>

Install Eigen::

    make install

.. note::
   Eigen is a header-only library, therefore there is nothing to "build" and install simply copies the header
   files. Installing via CMake however makes it easy for the library to be found automatically by EnDAS.
   If you omit the ``-DCMAKE_INSTALL_PREFIX=`` option, Eigen header files will be installed in the default 
   library location for your system.



Building the quick way
^^^^^^^^^^^^^^^^^^^^^^

Change (``cd``) to the EnDAS source directory and create subdirectory called ``build``::

    cd <endas-dir>
    mkdir build

Change (``cd``) to the build directory and run CMake::

    cd build
    cmake ..

.. note::
   EnDAS does not support in-source builds (i.e. without a build directory) and CMake will complain if run from
   the source directory.

This will configure the project to use default settings and the resulting libraries to be installed to the 
default install location for your system. Please refer to the :ref:`ENDAS_install_configure` section for the list of available 
configuration options.

At this point, CMake has generated all files needed to build EnDAS. On Linux, Makefiles are used by default. Staying in 
the ``build`` directory, run::

    make 

and wait for the library and examples to be built. If needed, the library can be installed on your system via::

    make install


Building and installing from source on Windows
----------------------------------------------

To be written...



.. _ENDAS_install_configure:

Configuring EnDAS
-----------------

EnDAS defines the following options in its CMake configuration:

- ``DISABLE_BLAS_LAPACK`` - Use ``-DDISABLE_BLAS_LAPACK=ON`` to disable use of BLAS and LAPACK
- ``DISABLE_OPENMP`` - Use ``-DDISABLE_OPENMP=ON`` to disable use of OpenMP
- ``DISABLE_PLOTTING`` - Use ``-DDISABLE_PLOTTING=ON`` to disable plotting functionality. This will
  remove the dependency on the Python 2.x development libraries
- ``DISABLE_PROFILING`` - Use ``-DDISABLE_PROFILING=ON`` to disable the integrated profiling functionality

























