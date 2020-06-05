/**
 * @file Variogram.hpp
 * @author Martin Gunia
 * 
 * Semivariograms, variograms and covariance functions.
 */

#ifndef __ENDAS_SPATIAL_VARIOGRAM_HPP__
#define __ENDAS_SPATIAL_VARIOGRAM_HPP__

#include <Endas/Config.h>
#include <Endas/Core/LinAlg.hpp>


namespace endas
{

/** 
 * @addtogroup spatial
 * @{ 
 */


/**
 * Covariance function of a stationary process.
 */
class ENDAS_DLL CovarianceFn
{
public:

    virtual ~CovarianceFn();

    /**
     * Returns values of the covariance function for two sets of locations ``A`` and ``B``.
     * The locations are stored column-wise so that the number of rows correponds to the number of
     * spatial dimensions and the number of columns corresponds to the number of locations.
     * The result is stored in the array ``out``.
     * 
     * @param A     d x n array of locations in set A where *d* is the number of spatial dimensions
     *              and *n* is the number of locations.
     * @param A     d x n array of locations in set B where *d* is the number of spatial dimensions
     *              and *n* is the number of locations.
     * @param out   Array where the result is stored.
     */
    virtual void values(const Ref<const Array2d> A, const Ref<const Array2d> B, 
                        Ref<Array> out) const = 0;
};


/**
 * Covariance function of a stationary isotropic process.
 * 
 * IsotropicCovarianceFn extends the CovarianceFn interface and subclasses only have to implement
 * the `IsotropicCovarianceFn::values(h, out)` method. The default implementation of 
 * `CovarianceFn::values(A, B, out)` assumes Euclidean coordinate system; if this is not desired, 
 * you can override this method as well.
 */
class ENDAS_DLL IsotropicCovarianceFn : public CovarianceFn
{
public:

    virtual ~IsotropicCovarianceFn();

    /**
     * Returns values of the covariance function at distances ``h``.
     * The result is stored in the array ``out``.
     * 
     * @param h     Array of distances for which the covariance value is returned
     * @param out   Array where the result is stored. The same array can be passed for ``r`` and 
     *              ``out``
     */
    virtual void values(const Ref<const Array> h, Ref<Array> out) const = 0;


    virtual void values(const Ref<const Array2d> A, const Ref<const Array2d> B, 
                        Ref<Array> out) const override;
};


/**
 * Covariance function from the exponential family. The covariance is defined by
 * 
 * @f[
 *   C(h) = \sigma e^{-(h/L)^{\alpha} }
 * @f]
 * 
 * where *h* is the distance and *L* is the correlation length.
 * 
 * Notable covariance functions of this family include the exponential covariance (``alpha=1``) 
 * and Gaussian covariance (``alpha=2``).
 */
class ENDAS_DLL ExponentialFamilyCovFn : public IsotropicCovarianceFn
{
public:

    /** 
     * ExponentialFamilyCovFn constructor.
     * 
     * @param alpha  The power term in the exponent.
     * @param L      Correlation length.
     * @param sigma  Standard deviation at r=0.
     * 
     */
    ExponentialFamilyCovFn(double alpha, double L, double sigma);

    virtual void values(const Ref<const Array> h, Ref<Array> out) const override;

private:
    double mAlpha;
    double mL;
    double mSigma;
};


/** 
 * Shorthand for ExponentialFamilyCovFn(1, sigma).
 * 
 * @param L      Correlation length.
 * @param sigma  Standard deviation at r=0.
 */
inline ExponentialFamilyCovFn ExponentialCovFn(double L, double sigma)
{
    return ExponentialFamilyCovFn(1, L, sigma);
}


/** 
 * Shorthand for ExponentialFamilyCovFn(2, sigma).
 * 
 * @param L      Correlation length.
 * @param sigma  Standard deviation at r=0.
 */
inline ExponentialFamilyCovFn GaussianCovFn(double L, double sigma)
{
    return ExponentialFamilyCovFn(2, L, sigma);
}


/**
 * Spherical covariance. The covariance is defined by the function
 * 
 * @f[
 *   C(r) =
       \begin{cases}
         \sigma \left( 1 - \left( \frac{2}{3}\frac{h}{L} - \frac{1}{2}\frac{h}{L}^3 \right) \right) & \mbox{if } h < L \\
        0  & \mbox{otherwise}
    \end{cases}
 * @f]
 * 
 * where *h* is the distance and *L* is the correlation length.
 */
class ENDAS_DLL SphericalCovFn : public IsotropicCovarianceFn
{
public:

    /** 
     * ExponentialFamilyCovFn constructor.
     * 
     * @param L      Correlation length.
     * @param sigma  Standard deviation at r=0.
     */
    SphericalCovFn(double L, double sigma);

    virtual void values(const Ref<const Array> h, Ref<Array> out) const override;

private:
    double mL;
    double mSigma;
};



/** @} */

}

#endif