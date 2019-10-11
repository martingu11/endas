Installation
============

Requirements
------------

The EnDAS library requires the following:

- `Python <https://www.python.org/>`_ 3.x (3.6 or newer recommended)
- `NumPy <https://numpy.org/>`_ and `SciPy <https://www.scipy.org/>`_
- Working C/C++ compiler for building the source distribution (see below)

Optionally, you may want to have the following installed as well:

- `Matplotlib <https://matplotlib.org/>`_ when running examples or plotting your own results
- `Sphinx <http://www.sphinx-doc.org>`_ to generate documentation

.. note::
   EnDAS is tested on CPython. It may be possible to use EnDAS with `PyPy<https://pypy.org/>`_ version 6.0.0 or above
   and SciPy 1.1.0 or newer. The instructions below assume you are using CPython.

Some parts of EnDAS are implemented as a C/C++ extension module (using Cython). Currently EnDAS comes in a
*source distribution* and you will therefore need a working C/C++ compiler in order to install it. On Linux, the
compiler is included with the OS and EnDAS should install out of the box. On Windows, you will need to install the
Microsoft Visual Studio (or the standalone compiler tools, if available). You should use Visual Studio version that
is recommended for your Python as listed `here <https://wiki.python.org/moin/WindowsCompilers>`_.


Installing from PyPI
--------------------

The easiest way to install EnDAS is to use the source distribution published in PyPi and install with ``pip``::

    pip install endas

This will build the C/C++ extension and install EnDAS into the Python environment form which you run ``pip``.


Building and installing from source
-----------------------------------

The source code is hosted on GitHub: `<https://github.com/martingu11/endas>`_. After cloning the repository, you can
install EnDAS in a working Python 3 distribution (with the dependencies above installed) by typing::

    python setup.py install --cython-rebuild

in the terminal from the root folder of the EnDAS library (where setup.py is located). The ``--cython-rebuild`` flag is
needed for the C/C++ extension source files to be generated from Cython sources. It needs to be done only once, unless
you are making changes to .pyx files.

If you plan on developing EnDAS or making changes to it, you may want to install it as and *ediable package*. This is
done with the "develop" command instead of "install"::

    python setup.py build_ext --inplace --cython-rebuild
    python setup.py develop

As above, this will build the extension modules and "install" EnDAS by creating a link to the source code. You do not
need to run ``setup.py develop`` afterwards. If you edit any Cython sources, you still need to rebuild the compiled
extension modules though::

    python setup.py build_ext --inplace --cython-rebuild


Building documentation
----------------------

If you have Sphinx installed (and the sphinx-build utility is on PATH), you can generate documentation by typing::

    cd doc
    make html

Please note that you also need GNU Make to run this (installed on all Linux flavours, for Windows version see
http://gnuwin32.sourceforge.net/packages/make.htm).


























