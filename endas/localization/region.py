"""
Region in N dimensions.

"""

__all__ = ['Region2d']


class Region2d:
    """
    Region in a two-dimensional XY plane defined by the intervals [x, xend) and [y, yend).
    """
    __slots__ = ('x', 'y', 'xend', 'yend')

    def __init__(self, x, y, xend, yend):
        self.x = x
        self.y = y
        self.xend = xend
        self.yend = yend

    def copy(self):
        return Region2d(self.x, self.y, self.xend, self.yend)

    def inflate(self, dx, dy):
        self.x -= dx
        self.y -= dy
        self.xend += dx
        self.yend += dy
        return self

    def clip(self, x, y, xend, yend):
        self.x = max(self.x, x)
        self.y = max(self.y, y)
        self.xend = min(self.xend, xend)
        self.yend = min(self.yend, yend)
        return self

    @property
    def shape(self): return self.yend - self.y, self.xend - self.x

    def area(self):
        shp = self.shape
        return shp[0] * shp[1]

    @property
    def centre(self):
        return self.x + (self.xend - self.x) / 2, self.y + (self.yend - self.y) / 2
