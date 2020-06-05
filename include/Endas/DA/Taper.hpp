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
 * @addtogroup da
 * @{ 
 */



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
 * The function is defined by
 *
 * @f[
 *   G(r) =
 *   \begin{cases}
 *     1 - \frac{5}{3}r^2 + \frac{5}{8}r^3 + \frac{1}{2}r^4 - \frac{1}{4}r^5 & \mbox{if } 0 \leq r < 1 \\
 *     4 - 5r + \frac{5}{3}r^2 + \frac{5}{8}r^3 - \frac{1}{2}r^4 + \frac{1}{12}r^5 - \frac{2}{3r} & \mbox{if } 0 \leq r < 2 \\
 *     0 & \mbox{if } r \geq 2
 *   \end{cases}
 * @f]
 *
 * where @f$ r = \vert i-j \vert / (L/2) @f$ and _i_, _j_ are the (indexes of) two state variables,
 * _L_ is the distance scale factor (the localization radius). Therefore, the support range of the
 * tapering function is _L_. 
 * 
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
 * Tapering function that does not perform any tapering at all.
 * 
 * The tapering function is defined as
 * 
 * @f[
 *   G(r) =
 *     \begin{cases}
 *       1 & \mbox{if } r < L \\
 *       0 & \mbox{if } r \geq L
 *     \end{cases}
 * @f]
 * 
 * The main purpose of NoTaper is to provide information about localization radius even if no
 * tapering is desired.
 */
class NoTaper : public TaperFn
{
public:

    /** 
     * Constructor.
     * 
     * @param L  The support range of the tapering function.
     */
    NoTaper(double L);
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
 * The taper function is defined by
 * 
 * @f[
 *   G(r) =
 *     \begin{cases}
 *       \left( 1 - \left( \frac{2}{3}\frac{r}{L} - \frac{1}{2}\frac{r}{L}^3 \right) \right) & \mbox{if } r < L \\
 *       0  & \mbox{if } r \geq L
 *     \end{cases}
 * @f]
 * 
 * where _r_ is the distance (between two state variables) and _L_ is the distance scale factor
 * (the localization radius). Therefore, the support range of the tapering function is _L_.
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


/** @} */

}


#endif