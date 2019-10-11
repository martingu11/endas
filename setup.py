import sys
import os.path
from setuptools import setup, Extension, find_packages

if sys.version_info[:2] < (3, 6):
    raise RuntimeError("Python >= 3.6 required.")

CYTHON_REBUILD = False
CYTHON_ANNOTATE = True


# Extra arguments to setup.py install/build/bukd_ext
# This is a bit hacky way to add arguments but avoids having to subclass the setuptools commands
# just to do something as simple as this.

# Force  re-compilation of .c files from .pyx with Cython
if "--cython-rebuild" in sys.argv:
    sys.argv.remove("--cython-rebuild")
    CYTHON_REBUILD = True

# Have Cython produce the HTML annotations when compiling .pyx files
if "--cython-annotate" in sys.argv:
    sys.argv.remove("--cython-annotate")
    CYTHON_ANNOTATE = True

# Locates all extension modules to compile and returns a list of setuptools.Extension instances.
# Our extensions are implemented in Cython, therefore there is one module per .c or .pyx file.
def find_extensions(dir, file_types, extlist=None):
    if extlist is None: extlist = []

    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        root, file_ext = os.path.splitext(path)
        file_ext = file_ext.lower()

        if os.path.isfile(path) and file_ext in file_types:
            ext_name = root.replace(os.path.sep, ".")
            extlist.append(Extension(ext_name, sources=[path]))
        elif os.path.isdir(path):
            extlist.extend(find_extensions(path, file_types, extlist))

    return extlist


# Building C/C++ extensions from the source distribution package
if not CYTHON_REBUILD:
    extensions = find_extensions('endas', ['.c'])

    if len(extensions) == 0:
        print("C/C++ source files not found. If building code from the source repository, please build with "
              "--cython-rebuild first to generate the source files.")
        sys.exit(1)
# Building using Cython from source .pyx files
else:
    try:
        from Cython.Build import cythonize
    except:
        print("Building with --cython-rebuild but Cython is not installed. Please install it first.")
        sys.exit(1)

    extensions = cythonize(find_extensions('endas', ['.pyx']), annotate=CYTHON_ANNOTATE)

with open('./README.md', 'r') as readme:
    long_description = readme.read()


setup(
    name="endas",
    version="0.1.0",
    description="Ensemble Data Assimilation library",
    long_description=long_description,
    packages=find_packages(exclude=['tests']), #  ["endas", "endas.algorithms", "endas.localization"],
    install_requires=['numpy', 'scipy'],
    ext_modules=extensions
)
