#ifndef __ENDAS_ALGORITHMS_ENKF_HPP__
#define __ENDAS_ALGORITHMS_ENKF_HPP__


#include <Endas/Caching/ArrayCache.hpp>
#include <Endas/Algorithm/Algorithm.hpp>

#include <memory>

namespace endas
{


/**
 * Base class for various Ensemble Kalman Filter/Smoother variants.
 * 
 * EnKFVariant implementations do most of the heavy-lifting of the EnsembleKalmanFilter
 * by computing the transform applied to the forecast ensemble applied during the analysis
 * update.
 */
class ENDAS_DLL EnKFVariant
{
    /**
     * Called before ensembleTransform() to compute additional data from the global ensemble.
     * 
     * 
     * The purpose of the method is to avoid passing the global ensemble to the ensembleTransform() 
     * method when domain localization is used. In most cases, the (local) observation operator is 
     * applied to the global ensemble, resulting in much smaller transformed ensemble which is passed
     * to ensembleTransform(), which may be executed in a separate process or on a different
     * machine altogether.
     * 
     * @param Ag    Global ensemble array (size nxN)
     * @param H     Observation operator instance. If analysis is localized, this is the localized 
     *              observation operator.
     */
    virtual void processElobalEnsemble(const Ref<const Array2d> Ag, const ObservationOperator& H) const;
    

    /** 
     * Computes analysis step as a transform applied to the ensemble.
     * 
     * Assuming the analysis update is computed as A_a = A_f * T, this gives the transform `T`. The 
     * result must be an NxN matrix, where N is the ensemble size. 
     */
    virtual void ensembleTransform(const Ref<const Array2d> A, const Ref<const Array> z, 
                                   const ObservationOperator& H, const CovarianceOperator& R, 
                                   double inflation, 
                                   Ref<Array2d> out) const = 0;

};

/**
 * Ensemble Kalman Filter and Smoother.
 * 
 * This is a generic Ensemble Kalman Filter/Smoother (EnKF/EnKS) implementation that covers many of 
 * the EnKF variants such as the stochastic EnKS and ensemble transform EnKS (ETKS). Both global and 
 * localized analysis is supported; the latter is performed by partitioning the state space into a 
 * set of disjoint local domains that are processed independently.
 */
class ENDAS_DLL EnsembleKalmanSmoother : public SequentialEnsembleSmoother
{
public:

    /** 
     * Constructs EnsembleKalmanSmoother.
     * 
     * @param variant   EnKF/EnKS variant to use. 
     * @param lag       Smoother lag (number of time steps).
     * @param cache     Array cache instance for storing intermediate filter solutions. If `nullptr` 
     *                  is given, arrays are cached in memory.
     * 
     * If `lag` is 0, only Kalman Filter solutions will be computed. Array caching is not used when 
     * only computing filter solutions. 
     */
    EnsembleKalmanSmoother(const EnKFVariant& variant, int lag = 0, std::shared_ptr<ArrayCache> cache = nullptr);

    ~EnsembleKalmanSmoother();

 
    virtual void beginSmoother(const Ref<const Array2d> E0, int k0) override;
    virtual void beginAnalysis(Ref<Array2d> E, int k) override;
    virtual void assimilate(const Ref<const Array> z, const ObservationOperator& H, 
                            const CovarianceOperator& R) override;
    virtual void endAnalysis() override;
    virtual void endSmoother() override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};



}


#endif