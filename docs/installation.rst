Installation
============

Requirements
------------

The EnDAS library requires the following:

- `Python <https://www.python.org/>`_ 3.x (3.6 or newer recommended)
- `NumPy <https://numpy.org/>`_ and `SciPy <https://www.scipy.org/>`_

Optionally, you may want to have the following installed as well

- `Matplotlib <https://matplotlib.org/>`_ when running examples or plotting your own results
- `Sphinx <http://www.sphinx-doc.org>`_ to generate documentation

Some parts of EnDAS are implemented as a C/C++ extension module (using Cython). Building of
the extension module is optional and you can skip it completely if you do not want to deal
with the additional complexity (EnDAS will fall back on a pure Python/NumPy implementation
then). The C/C++ implementation however offers better performance. At the moment the only way
to have the extension module installed is to build EnDAS from a source distribution and you will
need a working C/C++ compiler on your machine. The GCC compiler is typically included with all
Linux distributions. On Windows, you may want to install the Microsoft Visual Studio (Community
edition).


Installing from PyPI
--------------------

The easiest way to install EnDAS is to use the source distribution published in PyPi and
install with ``pip``::

    pip install endas

This will also build and install the C/C++ extension, therefore a working C/C++ compiler
is needed. If you do not want to build the C/C++ extension, tell ``pip`` to pass
``--no-cpp-ext`` to ``setup.py`` when installing EnDAS::

    pip install --install-option="--no-cpp-ext"


Building and installing from source
-----------------------------------

The source code is hosted on GitHub: `<https://github.com/martingu11/endas>`_. After
cloning the repository, you can install EnDAS in a working Python 3 distribution (with
the dependencies above installed) by typing::

    python setup.py install

in the terminal from the root folder of the EnDAS library (where setup.py is located).
You may need to run the command as a super user for the installation to succeed. As with
the PyPi installation, this will attempt to build the C/C++ extension before installing
the module. To disable this (for example if you do not have a working compiler), run::

    python setup.py install --no-cpp-ext

instead.

Building documentation
----------------------

If you have Sphinx installed (and the sphinx-build utility is on PATH), you can generate
documentation by typing::

    cd doc
    make html

Please note that you also need GNU Make to run this (installed on all Linux flavours,
for Windows version see http://gnuwin32.sourceforge.net/packages/make.htm).


























