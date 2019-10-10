"""
Domain localization on N-dimensional grids.
"""

__all__ = ['Grid2d']

import numpy as np

from . import SpatialQuery, TaperFn
from .bbox import BBox
from . import StateSpacePartitioning
from endas import _get_cython_impl


class Grid2d(StateSpacePartitioning):
    """
    Organizes state vector elements on a two-dimensional grid.

    The grid can be sparse, i.e. not all grid cells must be included in the state vector. Both cartesian and polar
    coordinates can be used, the grid itself is not tied to any particular coordinate system.

    Args:
        nx         : Number of cells in the grid along the *x* axis
        ny         : Number of cells in the grid along the *y* axis
        extent     : Extent of the grid given as :class:`BBox`. If ``None`` is passed, a region defined by
                     ``BBox((0, 0), (nx, ny))`` is assumed (i.e. the grid starts at (0, 0) and cell size is 1.
        cs         : Coordinate system of the grid. Must be an instance of :class:`endas.localization.CoordinateSystem`
                     and ``cs.ndim`` must be 2
        mask       : Flat array of indexes identifying grid cells included in the state vector. If ``None`` is passed
                     (this is the default), the grid is assumed to be dense.
        taper_fn   : Taper function that defines the localization radius. Must be an instance of
                     :class:`endas.localization.TaperFn` or ``None`` (see notes below).
        block_size : Size of each local domain along each dimension. The default value 1 results in analysis being
                     performed for each grid cell individually
        padding    : If greater than zero, each local domain will be padded by the given number of cells.

    If not all cells in the grid are a part of the state, the ``mask`` array must be given to define which grid cells
    correspond to which state vector elements. The mask array is a flat array of the same length as the state vector,
    with each element ``mask[i]`` containing the index of the cell corresponding to state variable :math:`x_i`. Cells
    in the grid are indexed from top to bottom in row-major order. If ``mask`` is ``None``, it is assumed that all cells
    are included in the state vector. How grid cells are mapped to the state vector via the mask array is shown below:

    .. figure:: ../figures/loc_grid2d_masking.png
       :scale: 90%

       Masking grid cells and their mapping to the state vector. Grid cells selected by the ``mask`` array
       are shown with gray background.

    .. note::
       It is valid to map several state variables to a single grid cell, for example if they contain different
       quantities in a multivariate setting. In this case the ``mask`` array must always be given, even if the grid is
       dense.

    """

    def __init__(self, nx : int, ny, extent, cs, mask=None, taper_fn=None, block_size=1, padding=0):
        self._nx = nx
        self._ny = ny
        self._mask = mask
        self._extent = extent
        self._cellsize = extent.size / (nx, ny)

        self._bs = block_size
        self._pad = padding
        self._domains = None


    def generate_domains(self):
        if self._domains is None:
            self._domains = []

            allcells = np.arange(self._ny*self._nx).reshape(self._ny, self._nx)
            full_bbox = BBox(min=(0, 0), end=(self._nx+1, self._ny+1), dtype=np.int)

            for y in range(0, self._ny, self._bs):
                for x in range(0, self._nx, self._bs):

                    d_box = BBox(min=(x, y), end=(x + self._bs, y + self._bs), dtype=np.int).intersect(full_bbox)
                    d_box_padded = d_box.copy()
                    d_box_padded = d_box_padded.inflate((self._pad, self._pad)).intersect(full_bbox)

                    # Select cells and corresponding variables from the state vector that are inside this block
                    d_cells = allcells[d_box_padded.y:d_box_padded.yend, d_box_padded.x:d_box_padded.xend]
                    d_svec = d_cells.ravel()[self._mask] if self._mask is not None else d_cells.ravel()

                    #Todo: Precompute blending mask for the padded region

                    if len(d_svec) > 0:
                        self._domains.append((d_box, d_svec))

        return self._domains


    def get_local_extent(self, domain):
        dbox = domain[0]
        return BBox(min=self._extent.min + dbox.min * self._cellsize,
                    end=self._extent.min + dbox.end * self._cellsize).intersect(self._extent)


    def get_local_state_size(self, domain):
        return len(domain[1])


    def get_local_state(self, domain, xg):
        raise NotImplementedError()


    def put_local_state(self, domain, xl, xg):
        raise NotImplementedError()




