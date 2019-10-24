cdef class BBox2d:
    cdef public double x, y, xend, yend

    cpdef double sizex(self)
    cpdef double sizey(self)
    cpdef BBox2d copy(self)

    cpdef void inflate(self, double dx, double dy)
    cpdef void intersect(self, BBox2d other)
    cpdef void clip(self, double x, double y, double endx, double endy)


cdef class BBox2i:
    cdef public int x, y, xend, yend

    cpdef int sizex(self)
    cpdef int sizey(self)
    cpdef BBox2i copy(self)

    cpdef void inflate(self, int dx, int dy)
    cpdef void intersect(self, BBox2i other)
    cpdef void clip(self, int x, int y, int endx, int endy)

