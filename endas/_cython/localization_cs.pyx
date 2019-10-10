"""
Cython implementation of endas/localization/cs.py
"""
import cython

from libc.math cimport sqrt, sin, cos, asin

ctypedef fused coord_t:
    int
    float
    double


# Euclidean distance specialized for the most common cases - 1d, 2d and 3d.
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_euclid_distance_1d(coord_t[:] A not None, coord_t[:] B not None, coord_t[:] out not None):
    # At this point we assume all the necessary asserts were already done! A and B are flat arrays
    cdef int n = A.shape[0]
    cdef int i
    cdef coord_t d
    for i in range(n):
        d = A[i] - B[i]
        if d < 0: d = -d   # Does anyone know better way to do abs() for both ints and floats?
        out[i] = d


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_euclid_distance_2d(double[:,:] A not None, double[:,:] B not None, double[:] out not None):
    # At this point we assume all the necessary asserts were already done! A and B have shape (n, 2)
    cdef int n = A.shape[0]
    cdef int i
    for i in range(n):
        out[i] = sqrt( (A[i,0]-B[i,0])**2 + (A[i,1]-B[i,1])**2 )


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_euclid_distance_3d(coord_t[:,:] A not None, coord_t[:,:] B not None, double[:] out not None):
    # At this point we assume all the necessary asserts were already done! A and B have shape (n, 3)
    cdef int n = A.shape[0]
    cdef int i
    for i in range(n):
        out[i] = sqrt( (A[i,0]-B[i,0])**2 + (A[i,1]-B[i,1])**2 + (A[i,2]-B[i,2])**2 )


# Euclidean distance for generic N-d case with extra inner loop
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_euclid_distance_Nd(coord_t[:,:] A not None, coord_t[:,:] B not None, double[:] out not None):
    # At this point we assume all the necessary asserts were already done! A and B have shape (n, N)
    cdef int n = A.shape[0]
    cdef int N = A.shape[1]
    cdef int i
    cdef double s
    for i in range(n):
        s = 0.0
        for j in range(N): s+= (A[i,j] - B[i,j])**2
        out[i] = sqrt(s)



# Lat-lon distance on a perfect sphere using the Haversine formula
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cs_latlon_distance(double[:,:] A not None, double[:,:] B not None, double[:] out not None, double R):
    # At this point we assume all the necessary asserts were already done! A and B have shape (n, 2)
    cdef int n = A.shape[0]
    cdef int i
    cdef double Alat, Blat, Alon, Blon
    cdef double dfif_lat, diff_lon
    cdef double a

    for i in range(n):
        Alat = A[i,0] * 0.0174532925    # Degree to Radian
        Alon = A[i,1] * 0.0174532925
        Blat = B[i,0] * 0.0174532925
        Blon = B[i,1] * 0.0174532925

        a = sin((Alat - Blat) / 2.0)**2 + cos(Alat) * cos(Blat) * sin((Alon - Blon) / 2.0)**2
        if a > 1.0: a = 1.0
        a = 2.0 * asin(sqrt(a))

        out[i] = a * R
