EnDAS internals implemented in Cython.

All Cython implementation of EnDAS routines should go here and be compiled into the
`endas._cython_impl` module (see setup.py). The compilation of the module is optional
on the user's side and Python/NumPy equivalents should be provided as well. EnDAS
code can call `endas._get_cython_impl()` to retrieve the `endas._cython_impl` module
(returns `None` if the module has not been compiled).
