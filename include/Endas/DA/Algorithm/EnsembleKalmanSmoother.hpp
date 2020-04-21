#ifndef __ENDAS_ALGORITHMS_ENSEMBLE_KALMAN_SMOOTHER_HPP__
#define __ENDAS_ALGORITHMS_ENSEMBLE_KALMAN_SMOOTHER_HPP__


#include <Endas/DA/Sequential.hpp>
#include <Endas/Caching/ArrayCache.hpp>

#include <memory>
#include <vector>

namespace endas
{


/**
 * Base class for various Ensemble Kalman Filter/Smoother variants.
 * 
 * EnKFVariant implementations do most of the heavy-lifting of the EnsembleKalmanFilter
 * by computing the transform applied to the forecast ensemble applied during the analysis
 * update.
 */
class ENDAS_DLL EnKSVariant
{
public:

    virtual ~EnKSVariant();

    /** Clones this instance. */
    virtual std::unique_ptr<EnKSVariant> clone() const = 0;

    /**
     * Called once before the variant is used in assimilation. 
     * This can be used to pre-compute data.
     */
    virtual void init(int n, int N);


    /**
     * Applies covariance inflation to the given ensemble.
     * 
     * @param E     Ensemble to apply inflation to.
     * @param f     Covariance inflation factor.
     * @param k     Time step index.
     * 
     * @note It is not necessary to inflate the ensemble `E`. Implementations may choose to 
     *       apply covariance inflation as part of the ensemble transform matrix. In this 
     *      case the coefficient `k` is typically just stored.
     * 
     * The default implementation does inflate the ensemble.
     */
    virtual void applyCovInflation(Ref<Array2d> E, double factor, int k) const;


    /**
     * Called before ensembleTransform() to compute additional data from the global ensemble.
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
     * 
     * @return Pointer to the data to be passed to ensembleTransform() (can be null pointer too).
     *
     * The default implementation returns unique_void_ptr(nullptr).
     */
    virtual void processGlobalEnsemble(const Ref<const Array2d> Eg, const ObservationOperator& H,
                                       int k, std::vector<Array2d>& dataOut) const;
    

    /** 
     * Computes analysis step as a transform applied to the ensemble.
     * 
     * Assuming the analysis update is computed as E_a = E_f * X, this gives the transform X. The 
     * ensemble `E` should be modified in-place and X should be stored in the pre-allocated Xout
     * array.
     */
    virtual void ensembleTransform(Ref<Array2d> E, std::vector<Array2d>& Egdata, 
                                   const Ref<const Array> z, const CovarianceOperator& R, int k,
                                   Ref<Matrix> Xout) const = 0;

};


/**
 * Ensemble Kalman Filter and Smoother.
 * 
 * This is a generic Ensemble Kalman Filter/Smoother (EnKF/EnKS) implementation that covers many 
 * EnKF/EnKS variants such as the stochastic Ensemble Kalman Smoother (EnKS) or variants based on the 
 * ensemble transform (ESTKS). See the documentation of the EnKSVariant class for information about 
 * implemented variants.
 * 
 * Both global and localized analysis is supported; the latter is performed by
 * partitioning the state space into a set of disjoint local domains that are processed independently.
 */
class ENDAS_DLL EnsembleKalmanSmoother : public SequentialEnsembleSmoother
{
public:

    /** 
     * Constructs EnsembleKalmanSmoother.
     * 
     * @param variant   EnKF/EnKS variant to use. The instance is copied.
     * @param lag       Smoother lag (number of time steps).
     * @param cache     Array cache instance for storing intermediate filter solutions. If `nullptr` 
     *                  is given, arrays are cached in memory.
     * 
     * If `lag` is 0, only Kalman Filter solutions will be computed. Array caching is not used when 
     * only computing filter solutions. 
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



//-------------------------------------------------------------------------------------------------
// Ensemble Kalman Filter/Smoother variants
//-------------------------------------------------------------------------------------------------


/**
 * Classic (stochastic) Ensemble Kalman Filter and Smoother variant with perturbed observations.
 */
class ENDAS_DLL EnKS : public EnKSVariant
{
public:

    virtual std::unique_ptr<EnKSVariant> clone() const override;

    virtual void processGlobalEnsemble(const Ref<const Array2d> Eg, const ObservationOperator& H,
                                       int k, std::vector<Array2d>& dataOut) const override;


    virtual void ensembleTransform(Ref<Array2d> E, std::vector<Array2d>& Egdata, 
                                   const Ref<const Array> z, const CovarianceOperator& R, int k, 
                                   Ref<Matrix> X5out) const override;
};



/**
 * Error Subspace Transform Kalman Filter and Smoother variant.
 */
class ENDAS_DLL ESTKS : public EnKSVariant
{
public:

    virtual std::unique_ptr<EnKSVariant> clone() const override;

    virtual void init(int n, int N) override;

    virtual void applyCovInflation(Ref<Array2d> E, double factor, int k) const override;

    virtual void processGlobalEnsemble(const Ref<const Array2d> Eg, const ObservationOperator& H,
                                       int k, std::vector<Array2d>& dataOut) const override;

    virtual void ensembleTransform(Ref<Array2d> E, std::vector<Array2d>& Egdata, 
                                   const Ref<const Array> z, const CovarianceOperator& R, int k, 
                                   Ref<Matrix> X5out) const override;

private:
    mutable double mInflation;
    Matrix mT;
};




}


#endif