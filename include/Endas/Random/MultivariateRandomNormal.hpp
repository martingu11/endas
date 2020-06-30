/**
 * @file MultivariateRandomNormal.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_RANDOM_MULTIVARIATE_RANDOM_NORMAL_HPP__
#define __ENDAS_RANDOM_MULTIVARIATE_RANDOM_NORMAL_HPP__

#include <Endas/Core/LinAlg.hpp>

#include <Eigen/LU>

namespace endas
{

/** 
 * @addtogroup random
 * @{ 
 */


/**
 * Random number generator sampling from a multivariate Normal distribution.
 * 
 * The generator is initialized with a covariance matrix and (optionally) mean vector.
 * To generate samples from the Normal distribution , use ``operator()``:
 * 
 *     MultivariateRandomNormal sampler(mean, cov);
 * 
 *     // Draw 100 samples from N(mean, cov)
 *     Array2d samples(mean.size(), 100);
 *     sampler(samples);
 * 
 * @rst
 * .. tip:: 
 *    If you do not need the covariance matrix or mean vector afterwards, you can pass them as 
 *    r-value references (using ``std::move()``) to re-use their memory for the multivariate 
 *    sampler:
 *        
 *    .. code-block:: cpp
 *  
 *       MultivariateRandomNormal sample(std::move(mean), std::move(cov));
 * 
 * @endrst
 * 
 */
class ENDAS_DLL MultivariateRandomNormal
{
public:

    /** 
     * MultivariateRandomNormal constructor, assuming zero mean.
     * 
     * @param cov   Covariance matrix of the mutivariate Normal distribution. 
     */
    MultivariateRandomNormal(Matrix cov);

    /** 
     * MultivariateRandomNormal constructor.
     * 
     * @param mean  Mean of the Normal distribution.
     * @param cov   Covariance matrix of the mutivariate Normal distribution.
     */
    MultivariateRandomNormal(Array mean, Matrix cov);


    ~MultivariateRandomNormal();

    /** 
     * Fills the given array with samples from the distribution.
     * 
     * If the ``out`` array has multiple columns, each column will contain an independent 
     * sample drawn from the distribution.
     */
    void sample(Ref<Array2d> out) const;


private:
    struct Data;
    std::unique_ptr<Data> mData;
   
};



/** @} */

}

#endif