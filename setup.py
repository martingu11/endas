import sys
import os.path
from setuptools import setup, Extension

if sys.version_info[:2] < (3, 6):
    raise RuntimeError("Python >= 3.6 required.")

CYTHON_REBUILD = False
CYTHON_ANNOTATE = True

# This is a bit hacky but will do the job and is the easiest way to add

# Force  recompilaton of .c files form .pyx with Cython
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


# Source files for the _cython_impl module containing Cython implementation of
# some parts of EnDAS.
cython_impl_sources = [
    "endas/_cython/localization_taper.pyx"
]

extensions = [
    # EnDAS internals implemented as Cython extension
    Extension("endas._cython_impl", finalize_sources(cython_impl_sources))
]


if CYTHON_REBUILD:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, annotate=CYTHON_ANNOTATE)



setup(
    name="endas",
    version="0.1",
    description="Ensemble Data Assimilation library",
    packages=["endas", "endas.algorithms", "endas.localization"],
    ext_modules=extensions
)
