"""
Domain localization on N-dimensional grids.
"""

__all__ = ['Grid2d']

import numpy as np
import math
import array

from . import SpatialQuery, TaperFn
from .bbox import BBox2i
from . import StateSpacePartitioning


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

    def __init__(self, nx, ny, extent, cs, mask=None, block_size=1, padding=0):
        self._nx = nx
        self._ny = ny
        self._mask = mask
        self._extent = extent
        self._cellsize = extent.size / (nx, ny)
        self._cs = cs

        self._bs = block_size
        self._pad = padding
        self._domains = None
        self._generate_domains()


    def _generate_domains(self):
        if self._domains is None:
            self._domains = []

            allcells = np.arange(self._ny*self._nx).reshape(self._ny, self._nx)
            full_bbox = BBox2i(0, 0, self._nx+1, self._ny+1)

            for y in range(0, self._ny, self._bs):
                for x in range(0, self._nx, self._bs):

                    d_box = BBox2i(x, y, x + self._bs, y + self._bs).intersect(full_bbox)
                    d_box_padded = d_box.copy()
                    d_box_padded.inflate(self._pad, self._pad)
                    d_box_padded.intersect(full_bbox)

                    # Select cells and corresponding variables from the state vector that are inside this block
                    d_cells = allcells[d_box_padded.y:d_box_padded.yend, d_box_padded.x:d_box_padded.xend]
                    d_svec = d_cells.ravel()[self._mask] if self._mask is not None else d_cells.ravel()

                    #Todo: Precompute blending mask for the padded region

                    if len(d_svec) > 0:
                        self._domains.append((d_box, d_svec))

        return self._domains


    def num_domains(self): return len(self._domains)


    def get_local_observations(self, domain_id, z_coords, taper_fn):

        assert domain_id >= 0 and domain_id <= self.num_domains
        d_box, _ = self._domains[domain_id]

        # Centre of the domain in "real-world" coordinates
        d_centre_coord = np.array((d_box.centerx() * self._cellsize + self._extent.x,
                                   d_box.centery() * self._cellsize + self._extent.y)).reshape(1,-1)

        r = int(math.ceil(taper_fn.support_range))

        # Grid2d allows the use of SpatialQuery or plain observation coordinates in `z_coords`. In the former
        # case the spatial index is used to find observations near the local domain. In the latter case a
        # brute-force search is done
        if isinstance(z_coords, np.ndarray):
            if z_coords.ndim != 2:
                raise ValueError("z_coords must be two-dimensional array")
            if z_coords.shape[1] != 2:
                raise ValueError("z_coords array must be of shape (n,2)")
            m = z_coords[0]

            dist = self._cs.distance(d_centre_coord, z_coords)
            selected = np.where(dist < r)[0]
            dist = dist[selected]
            return selected, dist

        # Assume z_coords is a SpatialQuery instance
        elif isinstance(z_coords, SpatialQuery):
            selected, dist = z_coords.range_query(d_centre_coord, r, distances=True)
            return selected, dist

        elif z_coords is None:
            raise ValueError("z_coords cannot be None")
        else:
            raise TypeError("z_coords is of unsupported type")


    def get_local_state_size(self, domain_id):
        assert domain_id >= 0 and domain_id <= self.num_domains
        return len(self._domains[domain_id][1])


    def get_local_state(self, domain_id, xg):
        assert domain_id >= 0 and domain_id <= self.num_domains
        raise NotImplementedError()


    def put_local_state(self, domain_id, xl, xg):
        assert domain_id >= 0 and domain_id <= self.num_domains
        raise NotImplementedError()






