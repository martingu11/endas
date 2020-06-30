#ifndef __ENDAS_ALGORITHMS_ENSEMBLE_KALMAN_SMOOTHER_HPP__
#define __ENDAS_ALGORITHMS_ENSEMBLE_KALMAN_SMOOTHER_HPP__


#include <Endas/DA/Sequential.hpp>
#include <Endas/DA/StateSpace.hpp>
#include <Endas/DA/StateSpacePartitioning.hpp>
#include <Endas/DA/Taper.hpp>
#include <Endas/Caching/ArrayCache.hpp>
#include "EnsembleKalmanSmootherVariant.hpp"

#include <memory>
#include <vector>


namespace endas
{

/** 
 * @addtogroup daalg
 * @{ 
 */


/**
 * Ensemble Kalman Filter and Smoother.
 * 
 * This is a generic Ensemble Kalman Filter/Smoother (EnKF/EnKS) implementation that covers many 
 * EnKF/EnKS variants. In fact, any EnKF/EnKS analysis scheme that can be expressed as a transform 
 * of the forecast ensemble @f$E_f@f$:
 * 
 * @f[
 *   E_a = E_f X
 * @f]
 * 
 * can be executed using this implementation. In addition, EnsembleKalmanSmoother also implements 
 * fixed-lag Kalman smoother via the sequential smoothing API (SequentialEnsembleSmoother). The update
 * step can operate on the entire state space (global analysis) or on local subsets of the state space
 * (i.e. localized analysis).
 */
class ENDAS_DLL EnsembleKalmanSmoother : public SequentialEnsembleSmoother
{
public:

    /** 
     * Constructs EnsembleKalmanSmoother.
     * 
     * If `lag` is 0, only Kalman Filter solutions will be computed. Array caching is not used when 
     * only computing filter solutions. 
     * 
     * @param variant   EnKF/EnKS variant to use. The instance is copied.
     * @param lag       Smoother lag (number of time steps).
     * @param cache     Array cache instance for storing intermediate filter solutions. If `nullptr` 
     *                  is given, arrays are cached in memory.
     */
    EnsembleKalmanSmoother(const EnKSVariant& variant, int n, int N, int lag = 0, 
                           std::shared_ptr<ArrayCache> cache = nullptr);

    EnsembleKalmanSmoother(const EnsembleKalmanSmoother&) = delete;
    ~EnsembleKalmanSmoother();

    /**
     * Sets the filter covariance inflation factor.
     */
    void setCovInflationFactor(double factor);


    /**
     * Sets the 'forgetting factor' of the smoother.
     */
    void setSmootherForgettingFactor(double factor);


    /**
     * Sets up analysis using localized domains.
     * 
     * The tapering function is used for reducing the influence of observations based on their 
     * distance from the local domain. If not given, the influence of observations is will not 
     * depend on the distance.
     * 
     * @param partitioner   State space partitioning implementation to use.
     * @param taperFn       Optional covariance tapering function.
     */
    virtual void localize(std::shared_ptr<const StateSpacePartitioning> partitioner, 
                          std::shared_ptr<const TaperFn> taperFn = nullptr);

    /**
     * Resets the smoother to a global analysis scheme.
     */
    virtual void globalize();

 
    virtual void beginSmoother(const Ref<const Array2d> E0, int k0) override;
    virtual void beginAnalysis(Ref<Array2d> E, int k) override;
    virtual void assimilate(const ObservationManager& omgr) override;
    virtual void endAnalysis() override;
    virtual void endSmoother() override;

private:
    
    struct Data;
    std::unique_ptr<Data> mData;
};



/** @} */

}


#endif