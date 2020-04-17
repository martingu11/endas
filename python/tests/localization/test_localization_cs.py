import pytest
import numpy as np
import time

from python.endas.localization import EuclideanCS, LatLonCS


def np_euclid_distance(a, b):
    diff = a - b
    diffSq = np.square(diff, out=diff)
    sqSum = np.sum(diffSq, axis=1)
    dist = np.sqrt(sqSum, out=sqSum)
    return dist


def np_latlon_distance(a, b, R):
    a *= 0.0174532925  # Degrees to radians
    b *= 0.0174532925
    A = np.sin((a[:,0] - b[:,0]) / 2.0) ** 2 + np.cos(a[:,0]) * np.cos(b[:,0]) * np.sin((a[:,1] - b[:,1]) / 2.0) ** 2
    A = np.clip(A, a_max=1.0, a_min=None, out=A)
    A = 2.0 * np.arcsin(np.sqrt(A))
    return A * R


@pytest.mark.parametrize("ndim", [1, 2, 3, 4, 5])
def test_euclid_cs(ndim):
    np.random.seed(1234)
    cs = EuclideanCS(ndim=ndim)

    print (ndim)

    # Try distance between 1000 random N-dimensional points, compared to pure NumPy implementation
    n = 1000
    a = np.random.randn(n, ndim)
    b = np.random.randn(n, ndim)

    start = time.process_time()
    dist_endas = cs.distance(a, b)
    end = time.process_time()
    print ("endas: {}ms".format((time.process_time() - start) * 1000.0))


    assert dist_endas.ndim == 1
    assert dist_endas.shape == (n,)

    start = time.process_time()
    dist_np = np_euclid_distance(a, b)
    end = time.process_time()
    print ("numpy: {}ms".format((time.process_time() - start) * 1000.0))

    assert np.allclose(dist_endas, dist_np)

    # Try with `a` having only 1 row. This is allowed by EnDAS and we should get distances from points
    # in `b` from the single point in `a`
    a1 = np.random.randn(1, ndim)
    dist_endas = cs.distance(a1, b)
    assert dist_endas.ndim == 1
    assert dist_endas.shape == (n,)

    dist_np = np_euclid_distance(a1, b)
    assert np.allclose(dist_endas, dist_np)



def test_latlon_cs():
    np.random.seed(1234)
    cs = LatLonCS()

    # Try distance between 1000 random points, compared to pure NumPy implementation
    n = 1000
    a = np.hstack((
        np.random.uniform(-90.0, 90.0, n).reshape(n, 1),
        np.random.uniform(-180.0, 180.0, n).reshape(n, 1)
    ))
    b = np.hstack((
        np.random.uniform(-90.0, 90.0, n).reshape(n, 1),
        np.random.uniform(-180.0, 180.0, n).reshape(n, 1)
    ))

    dist_endas = cs.distance(a, b)
    assert dist_endas.ndim == 1
    assert dist_endas.shape == (n,)

    dist_np = np_latlon_distance(a, b, cs.R)
    assert np.allclose(dist_endas, dist_np)

    # Try with `a` having only 1 row. This is allowed by EnDAS and we should get distances from points
    # in `b` from the single point in `a`
    a1 = np.hstack((
        np.random.uniform(-90.0, 90.0, 1),
        np.random.uniform(-180.0, 180.0, 1)
    )).reshape(1, -1)

    dist_endas = cs.distance(a1, b)
    assert dist_endas.ndim == 1
    assert dist_endas.shape == (n,)

    dist_np = np_latlon_distance(a1, b, cs.R)
    assert np.allclose(dist_endas, dist_np)







