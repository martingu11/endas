"""
Domain localization on N-dimensional grids.



"""

__all__ = ['GridLocalization2d']

import numpy as np

from endas import SpatialQuery
from endas.localization.region import Region2d
from endas.localization.dl_base import DomainLocalizationBase
from endas import _get_cython_impl


class GridLocalization2d(DomainLocalizationBase):
    """
    Implements domain localization of the analysis update on a 2-dimensional (sparse) grid.

    This localization implementation assumes that the state vector variables are formed from a two-dimensional grid.
    The grid can be fully or partially represented, i.e. not every grid cell must be a part of the state vector.


    Args:
        cs : The coordinate system of the grid. Must be an instance of :class:`endas.CoordinateSystem` and
             ``cs.ndim`` must be equal to 2

    """

    def __init__(self, cs, nx, ny, sv_query, extent=None, blocksize=3):
        super().__init__(cs)
        self._nx = nx
        self._ny = ny
        self._extent = extent if extent is not None else Region2d(0, 0, nx, ny)

        assert isinstance(sv_query, SpatialQuery)
        self._sv_query = sv_query

        self._bs = blocksize
        self._pad = blocksize
        self._domains = None

    def generate_domains(self):
        if self._domains is None:
            self._domains = []

        for y in range(0, self._ny, self._bs):
            for x in range(0, self._nx, self._bs):

                rect = Region2d(x, y, x + self._bs, y + self._bs).clip(0, 0, self._nx, self._ny)
                padded = rect.copy()
                padded = padded.inflate(self._pad, self._pad).clip(0, 0, self._nx, self._ny)

                # Select states that are actually within the window
                self._sv_query.range_query(lbound=(padded.x, padded.y), (padded.xend-1, padded.yend-1))



                if self._scoord is not None:
                    dist = dist2d_from_interval(self._scoord[:, 0], self._scoord[:, 1],
                                                rect.x, rect.xend - 1,
                                                rect.y, rect.yend - 1,
                                                self._pad, self._pad ** 2,
                                                out=dist)

                    selstates = np.where(dist <= self._pad)[0]
                    # Nothing in domain -> skip
                    if len(selstates) == 0: continue

                    seldist = dist[selstates]
                    weights = 1.0 - (seldist / self._pad) if self._pad > 0 else None

                else:
                    selstates = None
                    # Construct blending mask for the domain. If mask for this tile shape does not exist,
                    # it will be created and cached
                    # bm = blendmasks[rect.shape]
                    # bmx = rect.x - padded.x
                    # bmy = rect.y - padded.y
                    # bmxend = bm.shape[1] - (padded.xend - rect.xend)
                    # bmyend = bm.shape[0] - (padded.yend - rect.yend)
                    # bm = bm[bmy:bmyend, bmx:bmxend]

                self._domains.append((rect, padded, selstates, weights))

        return self._domains



    def get_state(self, domain, xg):
        raise NotImplementedError()

    def put_state(self, domain, xl, xg):
        raise NotImplementedError()




