#ifndef __ENDAS_ALGORITHMS_ENSEMBLE_KALMAN_SMOOTHER_VARIANT_HPP__
#define __ENDAS_ALGORITHMS_ENSEMBLE_KALMAN_SMOOTHER_VARIANT_HPP__

#include <Endas/DA/ObservationOperator.hpp>
#include <Endas/DA/CovarianceOperator.hpp>

#include <memory>
#include <vector>

namespace endas
{


/**
 * Base class for various Ensemble Kalman Filter/Smoother variants.
 * 
 * EnKFVariant implementations do most of the heavy-lifting of the EnsembleKalmanFilter by
 * computing the transform applied to the forecast ensemble applied during the analysis
 * update.
 * 
 * @rst
 * .. note::
 *      The instance is cloned for each executing thread and it is therefore not necessary 
 *      to implement locking.
 * @endrst
 */
class ENDAS_DLL EnKSVariant
{
public:

    virtual ~EnKSVariant();

    /** Clones this instance. */
    virtual std::unique_ptr<EnKSVariant> clone() const = 0;

    /**
     * Called once before the instance is used in assimilation. 
     * This can be used to pre-compute data.
     * 
     * @param n     State size
     * @param N     Ensemble size    
     */
    virtual void init(int n, int N);


    /**
     * Applies covariance inflation to the given ensemble.
     * 
     * The default implementation inflates the ensemble by the given factor. If covariance inflation 
     * can included exactly in the exsemble transform, the inflation factor ``f`` can be stored and 
     * used in ensembleTransform().
     * 
     * @param E     Ensemble to apply inflation to.
     * @param f     Covariance inflation factor.
     * @param k         The analysis time step index
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
     * @param k     The analysis time step index 
     * 
     * @return Pointer to the data to be passed to ensembleTransform() (can be null pointer too).
     */
    virtual void processGlobalEnsemble(const Ref<const Array2d> Eg, const ObservationOperator& H,
                                       int k, std::vector<Array2d>& dataOut) const;
    

    /** 
     * Computes analysis step as a transform applied to the ensemble.
     * 
     * This computes the analysis update as a transform of the forecast ensemble @f$ E_a = E_f * X @f$.
     * The ensemble is modified in-place and the ensemble transform *X* is stored in ``Xout``.
     * 
     * @param E         Ensemble to transform (nxN array)
     * @param Egdata    Data computed from the global ensemble by processGlobalEnsemble()    
     * @param z         Array of observed values
     * @param k         The analysis time step index
     * @param Xout      Pre-allocated matrix of size NxN where the ensemble transform is stored
     */
    virtual void ensembleTransform(Ref<Array2d> E, std::vector<Array2d>& Egdata, 
                                   const Ref<const Array> z, const CovarianceOperator& R, int k,
                                   Ref<Matrix> Xout) const = 0;

};


//-------------------------------------------------------------------------------------------------
// Ensemble Kalman Filter/Smoother variants
//-------------------------------------------------------------------------------------------------

/**
 * Classic (stochastic) Ensemble Kalman Filter/Smoother.
 * 
 * EnKS implements traditional Ensemble Filter/Smoother with perturbed observations. 
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
 * Error Subspace Transform Kalman Filter/Smoother.
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