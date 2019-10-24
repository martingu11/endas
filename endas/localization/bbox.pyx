"""
Axis-aligned bounding boxes in various dimensions.
"""

__all__ = [ 'BBox2d', 'BBox2i']

#import cython
#cimport cython


cdef class BBox2d:
    """
    Axis-aligned two-dimensional bounding box.

    The bounding box is defined by the half-open intervals [`x``, ``xend``) and [`y``, ``yend``), so the lower bound
    is inclusive while the upper bound is exclusive. This is to make it possible to define two touching bounding boxes
    in :math:`\\mathbf{R}^2` that are nevertheless disjoint.
    """

    def __cinit__(self, double x, double y, double xend, double yend):
        self.x = x
        self.y = y
        self.xend = xend
        self.yend = yend

    @property
    def shape(self):
        return self.xend - self.x, self.yend - self.y

    cpdef double sizex(self):
        return self.xend - self.x

    cpdef double sizey(self):
        return self.yend - self.y

    cpdef BBox2d copy(self):
        """
        Deep copy of the bounding box.
        """
        cdef BBox2d bb = BBox2d.__new__(BBox2d)
        bb.x = self.x
        bb.y = self.y
        bb.xend = self.xend
        bb.yend = self.yend
        return bb

    cpdef void inflate(self, double dx, double dy):
        self.x -= dx
        self.y -= dy
        self.xend += dx
        self.yend += dy

    cpdef void intersect(self, BBox2d other):
        self.clip(other.x, other.y, other.xend, other.yend)

    cpdef void clip(self, double x, double y, double xend, double yend):
        self.x = max(self.x, x)
        self.y = max(self.y, y)
        self.xend = min(self.xend, xend)
        self.yend = min(self.yend, yend)



cdef class BBox2i:
    """
    Axis-aligned two-dimensional bounding box.

    The bounding box is defined by the half-open intervals [`x``, ``xend``) and [`y``, ``yend``), so the lower bound
    is inclusive while the upper bound is exclusive. This is to make it possible to define two touching bounding boxes
    in :math:`\\mathbf{R}^2` that are nevertheless disjoint.
    """

    def __cinit__(self, int x, int y, int xend, int yend):
        self.x = x
        self.y = y
        self.xend = xend
        self.yend = yend

    @property
    def shape(self):
        return self.xend - self.x, self.yend - self.y

    cpdef int sizex(self):
        return self.xend - self.x

    cpdef int sizey(self):
        return self.yend - self.y

    cpdef BBox2i copy(self):
        """
        Deep copy of the bounding box.
        """
        cdef BBox2i bb = BBox2i.__new__(BBox2i)
        bb.x = self.x
        bb.y = self.y
        bb.xend = self.xend
        bb.yend = self.yend
        return bb

    cpdef void inflate(self, int dx, int dy):
        self.x -= dx
        self.y -= dy
        self.xend += dx
        self.yend += dy

    cpdef void intersect(self, BBox2i other):
        self.clip(other.x, other.y, other.xend, other.yend)

    cpdef void clip(self, int x, int y, int xend, int yend):
        self.x = max(self.x, x)
        self.y = max(self.y, y)
        self.xend = min(self.xend, xend)
        self.yend = min(self.yend, yend)




