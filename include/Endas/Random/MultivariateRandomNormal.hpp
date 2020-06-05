/**
 * @file MultivariateRandomNormal.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_RANDOM_MULTIVARIATE_RANDOM_NORMAL_HPP__
#define __ENDAS_RANDOM_MULTIVARIATE_RANDOM_NORMAL_HPP__

#include <Endas/Core/LinAlg.hpp>

namespace endas
{

/** 
 * @addtogroup random
 * @{ 
 */


/**
 * Random number generator sampling from a multivariate Normal distribution.
 */
class ENDAS_DLL MultivariateRandomNormal
{
public:

    /** 
     * MultivariateRandomNormal constructor, assuming zero mean.
     * 
     * @param cov   Covariance matrix of the mutivariate Normal distribution.
     */
    MultivariateRandomNormal(const Ref<const Matrix> cov);

    /** 
     * MultivariateRandomNormal constructor.
     * 
     * @param mean  Mean of the Normal distribution.
     * @param cov   Covariance matrix of the mutivariate Normal distribution.
     */
    MultivariateRandomNormal(const Ref<const Array> mean, const Ref<const Matrix> cov);

    /** 
     * Fills the given array with samples from the distribution.
     */
    void operator()(Ref<Array2d> out) const;



private:
    
    mutable Array mX;
    Array mMean;
    Matrix mTransform;
};



/** @} */

}

#endif