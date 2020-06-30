/**
 * @file CoordinateSystem.hpp
 * @author Martin Gunia
 * 
 * Coordinate systems.
 */

#ifndef __ENDAS_SPATIAL_COORDINATE_SYSTEM_HPP__
#define __ENDAS_SPATIAL_COORDINATE_SYSTEM_HPP__

#include <Endas/Core/LinAlg.hpp>

namespace endas
{

/** 
 * @addtogroup spatial
 * @{ 
 */


/**
 * Coordinate system abstraction.
 * 
 * Coordinate system is characterized by its dimension (i.e. the number of components each
 * coordinate takes) and must provide a distance metric via the `distance()` method.
 */ 
class ENDAS_DLL CoordinateSystem
{
public:

    virtual ~CoordinateSystem();


    /**
     * Returns the number of spatial dimensions.
     */
    virtual int dim() const = 0;

    /**
     * Returns `true` if this is a Cartesian coordinate system.
     */
    virtual bool isCartesian() const = 0;

    /**
     * Returns the distance between pairs of points.
     * 
     * The method computes distance between pairs of points in the sets `A` and `B`. The points are stored 
     * column-wise and their dimension must equal `ndim()` (therefore both `A` and `B` must be `ndim()` x `n` 
     * arrays where `n` is the number of points to process). Alternatively, the `B` array can contain a single
     * column, in this case the distances are calculated between all points in `A` and the single point in `B`.
     *
     * @param A     First set of points.
     * @param B     Second set of points.
     * @param out   Array of size `n` where to store the distances.
     */
    virtual void distance(const Ref<const Array2d> A, const Ref<const Array2d> B,
                          Ref<Array> out) const = 0;

};



/** 
 * Euclidean coordinate system in N dimensions.
 */
class ENDAS_DLL EuclideanCS : public CoordinateSystem
{
public:

    /**
     * EuclideanCS constructor.
     * 
     * @param dim  Number of dimensions, must be >= 1.
     */
    EuclideanCS(int dim);

    virtual int dim() const override;
    virtual bool isCartesian() const override;
    virtual void distance(const Ref<const Array2d> A, const Ref<const Array2d> B,
                          Ref<Array> out) const override;

private:
    int mNdim;
};



/** 
 * Polar coordinate system on a perfect sphere.
 * 
 * LatLonCS implements coordinate system on a perfect sphere with coordinates of any point expressed 
 * as latitude and longitude. The currently implemented coordinate system is strictly two-dimensional, 
 * i.e. there is no vertical component. The coordinates are assumed to be in degrees and the first 
 * coordinate is the latitude, second the longitude.
 * 
 * This class is a simple implementation of a polar coordinate system assuming a perfect sphere of 
 * radius `R` and uses the `Haversine formula <https://en.wikipedia.org/wiki/Haversine_formula>`_ to 
 * compute the great-circle distance. If the use of a better approximation is required (i.e. a 
 * spheroid), please consider coding your own.
 */
class ENDAS_DLL LatLonCS : public CoordinateSystem
{
public:

    /**
     * LatLonCS constructor.
     * 
     * @param R  The radius of the great circle on which the distance is calculated. The default value 
     *           of 6371 km corresponds to the mean Earth radius (R1) as defined by the International 
     *           Union of Geodesy and Geophysics.
     */
    LatLonCS(double R = 6.371e6);

    virtual int dim() const override;
    virtual bool isCartesian() const override;
    virtual void distance(const Ref<const Array2d> A, const Ref<const Array2d> B,
                          Ref<Array> out) const override;

private:
    int mR;
};


/** @} */

}

#endif