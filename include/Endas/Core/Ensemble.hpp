/**
 * @file Ensemble.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_DA_ENSEMBLE_HPP__
#define __ENDAS_DA_ENSEMBLE_HPP__

#include <stddef.h>

#include <Endas/Core/LinAlg.hpp>
#include <Endas/DA/CovarianceOperator.hpp>


namespace endas
{

/** 
 * @addtogroup core
 * @{ 
 */


/** 
 * Generates new ensemble with given mean and covariance.
 * 
 * @param u    Ensemble mean.
 * @param cov  Ensemble covariance.
 * @param out  Array where to store the result. Must be of size `u.size() x N, where N is 
 *             the desired number of ensemble members to generate.
 */
ENDAS_DLL void generateEnsemble(const Ref<const Array> u, const CovarianceOperator& cov, Ref<Array2d> out);


/**
 * Returns ensemble mean as an Eigen expression. 
 */ 
template <class DerivedIn>
inline auto ensembleMean(const Eigen::ArrayBase<DerivedIn>& E) -> decltype ( E.rowwise().mean() )
{
    return E.rowwise().mean();
} 



/** 
 * Transforms ensemble state vectors to anomaly (deviatiom from ensemble mean).
 * 
 * @param E    Ensemble array to transform.
 * @param out  Array where to store the result. Can be same as `E` (transformed in-place).
 */
ENDAS_DLL void toAnomaly(const Ref<const Array2d> E, Ref<Array2d> out);



/** 
 * Inflates an ensemble around its mean by given factor.
 * 
 * @param E  Ensemble array to inflate.
 * @param k  The inflation factor.
 */
ENDAS_DLL void inflateInPlace(Ref<Array2d> E, double k);


/** 
 * Inflates an ensemble around its mean by given factor.
 * This overload avoids calculaiton of the ensemble mean if it is already known.
 * 
 * @param E  Ensemble array to inflate.
 * @param k  The inflation factor.
 * @param Eu The ensemble mean. 
 */
ENDAS_DLL void inflateInPlace(Ref<Array2d> E, double factor, Ref<Array> Eu);



/** @} */

}

#endif