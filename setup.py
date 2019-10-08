import sys
import os.path
from setuptools import setup, Extension

if sys.version_info[:2] < (3, 6):
    raise RuntimeError("Python >= 3.6 required.")

BUILD_CPP_EXT = True
CYTHON_REBUILD = False
CYTHON_ANNOTATE = True


# Extra arguments to setup.py install/build/bukd_ext
# This is a bit hacky way to add arguments but avoids having to subclass the setuptools commands
# just to do something as simple as this.

# Disable building of the C/C++ extension
if "--no-cpp-ext" in sys.argv:
    sys.argv.remove("--no-cpp-ext")
    BUILD_CPP_EXT = True

# Force  recompilaton of .c files from .pyx with Cython
if "--cython-rebuild" in sys.argv:
    sys.argv.remove("--cython-rebuild")
    CYTHON_REBUILD = True

# Have Cython produce the HTML annotations when compiling .pyx files
if "--cython-annotate" in sys.argv:
    sys.argv.remove("--cython-annotate")
    CYTHON_ANNOTATE = True



def finalize_sources(sources):
    if not CYTHON_REBUILD:
        return [ os.path.splitext(s)[0] + ".c" for s in sources ]
    else:
        return sources


# Source files for the _cython_impl module containing Cython implementation of some parts of EnDAS.
cython_impl_sources = [
    "endas/_cython/localization_taper.pyx",
    "endas/_cython/localization_cs.pyx"
]

extensions = []

# Building EnDAS internals implemented as C/C++ extension
if BUILD_CPP_EXT:
    extensions.append(Extension("endas.cython_impl",
                                finalize_sources(cython_impl_sources)
                                ))


if BUILD_CPP_EXT and CYTHON_REBUILD:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, annotate=CYTHON_ANNOTATE)



setup(
    name="endas",
    version="0.1",
    description="Ensemble Data Assimilation library",
    packages=["endas", "endas.algorithms", "endas.localization"],
    install_requires=['numpy', 'scipy'],
    ext_modules=extensions,

)
