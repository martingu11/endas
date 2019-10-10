"""
Bounding boxes.
"""

__all__ = [ 'BBox']

import numpy as np


class BBox:
    """
    Axis-aligned N-dimensional bounding box.

    The bounding box is defined by the half-open interval [``min``, ``end``), so the lower bound is inclusive while
    the upper bound is exclusive. This is to make it possible to define two touching bounding boxes in
    :math:`\\mathbf{R}^N` that are nevertheless disjoint.

    """
    __slots__ = ('_N', '_min', '_end')   # _end rather than _max since it is open

    def __init__(self, ndim=None, min=None, end=None, dtype=np.double):
        if min is not None or end is not None:
            assert min is not None and end is not None
            assert len(min) == len(end)

        # Todo: This is not terribly efficient!
        self._N = ndim if ndim is not None else len(min)
        self._min = np.zeros(self._N, dtype=dtype) if dtype is not None else None
        self._end = np.zeros(self._N, dtype=dtype) if dtype is not None else None
        if min is not None: self._min[:] = min
        if end is not None: self._end[:] = end


    @property
    def ndim(self):
        """
        Number of dimensions of the bounding box.
        """
        return self._N

    @property
    def min(self): return self._min

    @property
    def end(self): return self._end

    def copy(self):
        """
        Deep copy of the bounding box.
        """
        bb = BBox(self._N, dtype=None)
        bb._min = np.copy(self._min)
        bb._end = np.copy(self._end)

    def inflate(self, d):
        self._min-= d
        self._end+= d
        return self


    def intersect(self, other):
        self.clip(other.min, other.end)
        return self


    def clip(self, min, end):
        self._min = np.clip(self._min, a_min=min, a_max=None, out=self._min)
        self._end = np.clip(self._end, a_min=None, a_max=end, out=self._end)
        return self

    @property
    def size(self): return self._end - self._min


