/**
 * @file Taper.hpp
 * @author Martin Gunia
 * 
 * Covariance tapering functions.
 */

#ifndef __ENDAS_DA_TAPER_HPP__
#define __ENDAS_DA_TAPER_HPP__

#include <Endas/Core/Core.hpp>


namespace endas
{


/**
 * Covariance tapering function with local support.
 * 
 * Tapering functions are used for localization of the analysis update by adjusting the influence
 * of observations based on their distance from the local state variable.
 */ 
class TaperFn
{
public:


    TaperFn(double L);

    virtual ~TaperFn();

    /**
     * Returns the support range of the tapering function, i.e. the distance `d` at which the 
     * tapering coefficient becomes zero.
     */
    virtual double supportRange() const;

    /**
     * Tapers the vector `x` based on corresponding distances `d`.
     * 
     * Tapering is done by multiplying each element `x_i` in the array `x` by weight `w(d_i)`,
     * where `d_i` is the distance assigned to the element `x_i`.
     * 
     * @param x     Array of values to taper
     * @param d     Array of distances 
     * @param out   Array where to store the result
     * 
     * The `x` array can be passed via `out` to perform the tapering in-place. 
     */
    virtual void taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const = 0;

protected:
    double mL;
};


/**
 * Gaspari-Cohn covariance tapering function.
 * 
 * The support range of the function is `L`. 
 */
class GaspariCohnTaper : public TaperFn
{
public:

    /** 
     * Constructor.
     * 
     * @param L  The support range of the tapering function.
     */
    GaspariCohnTaper(double L);
    virtual void taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const override;
};


/**
 * Linear covariance tapering function.
 * 
 * The tapering function implements linear falloff towards zero at its support range `L`.
 */
class LinearTaper : public TaperFn
{
public:

    /** 
     * Constructor.
     * 
     * @param L  The support range of the tapering function.
     */
    LinearTaper(double L);
    virtual void taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const override;
};


/**
 * Spherical covariance tapering function.
 * 
 */
class SphericalTaper : public TaperFn
{
public:

    /** 
     * Constructor.
     * 
     * @param L  The support range of the tapering function.
     */
    SphericalTaper(double L);
    virtual void taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const override;
};




}

#endif