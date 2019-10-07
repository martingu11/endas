# EnDAS Installation

 ## Requirements

The EnDAS library requires the following:

- Python 3.x (3.6+ recommended)
- NumPy and SciPy Python packages

Optionally, you may want to have the follwing installed as well

- Matplotlib Python package
- Sphinx (http://www.sphinx-doc.org/en/master) to generate documentation

To build and install EnDAS from source, you will also need:

- Cython


## Installing from PyPI

TO BE COMPLETED!


## Building and installing from source 

You can install EnDAS in a working Python 3 distribution (with the dependencies above installed)
by typing

    python3 setup.py install
   
in the terminal from the root folder of the EnDAS library (where `setup.py` is located). 
You may need to run the command as a super user for the installation to succeed. 


If you have Sphinx installed (and the `sphinx-build` utility is on `PATH`), you can generate 
documentation by typing

    cd doc
    make html

Please note that you also need GNU Make to run this (installed on all Linux flavours, 
for Windows version see http://gnuwin32.sourceforge.net/packages/make.htm).  


## Usage

Examples are provided in the ``examples`` directory. Each example file provides instructions
on how it should be run although typically it is as simple as just executing the Python file.



