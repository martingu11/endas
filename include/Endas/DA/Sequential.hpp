/**
 * @file Sequential.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_ALGORITHM_ALGORITHM_HPP__
#define __ENDAS_ALGORITHM_ALGORITHM_HPP__

#include "Model.hpp"
#include "CovarianceOperator.hpp"
#include "ObservationOperator.hpp"

#include <functional>
#include <limits>


namespace endas
{

/**
 * Special lag value indicating that a fixed-interval Kalman Smoother should be used, if available.
 * If not supported, the Kalman Smoother implementation will fall back on fixed-lag implementation
 * with a very large lag value.
 */
constexpr int LAG_FIKS = std::numeric_limits<int>::max();


/**
 * Sequential filter interface.
 */
class ENDAS_DLL SequentialFilter
{
public:

    /** 
     * Function called from endAnalysis() when solution is available. 
     */
    typedef std::function<void(const Ref<const Array> x, const Ref<const Matrix> P, int k)> OnResultFn;

    SequentialFilter();
    ~SequentialFilter();
    SequentialFilter(const SequentialFilter&) = delete;


    /**
     * Sets the callable that will be called every time a filter or smoother solution is available.
     * 
     * If lag is 0, `fn` is called on the filter solutions. More usefully, when lag > 0, `fn` will be 
     * called for every smoother solution that becomes available.
     * 
     * @param fn    Any callable compatible with OnResultFn.
     * 
     * The callable is copied.
     */
    virtual void onResult(OnResultFn fn);


    /**
     * Implements the forecast step.
     * 
     * The state vector `x` and covariance `P` are modified in-place and the system is propagated from 
     * time `t` to `t+dt`.
     * 
     * @param x     The state vector.
     * @param P     State error covariance matrix.
     * @param Q     Model error covariance matrix. Can be empty matrix for a perfect model.
     * @param k     The time stepping index passed to the evolution model.
     * @param dt    The time stepping interval.
     */
    virtual void forecast(Ref<Array> x, Ref<Matrix> P, const Ref<const Matrix> Q, int k, double dt) = 0;


    /**
     * Initiates the analysis/update step.
     * 
     * The reference to the state vector `x` is stored and will be updated as observations are being 
     * assimilated. 
     * 
     * @param x     The state vector at time `t`.
     * @param P     State error covariance matrix before assimilation of observations.
     * @param k     The time stepping index.
     * 
     */ 
    virtual void beginAnalysis(Ref<Array> x, Ref<Matrix> P, int k) = 0;

    /**
     * Assimilates observations in the vector `z` to the state vector `x`.
     * 
     * @param z     Array of observations.
     * @param H     Observation operator as a matrix.
     * @param R     Observation error covariance matrix.
     */
    virtual void assimilate(const Ref<const Array> z, const Ref<const Matrix> H, const Ref<const Matrix> R) = 0;
    
    /**
     * Must be called after all observations have been assimilated.
     */ 
    virtual void endAnalysis() = 0;


protected:
    OnResultFn mOnResultFn;

};


/**
 * Sequential smoother interface.
 */
class ENDAS_DLL SequentialSmoother : public SequentialFilter
{
public:

    /** 
     * Must be called first to provide the initial state vector and error covariance.
     * 
     * @param x0    Initial state vector at time step `k==0`
     * @param P0    Initial error covariance at time step `k==0`   
     * @param k0    The initial time step index (typically 0)
     */
    virtual void beginSmoother(const Ref<const Array> x0, const Ref<const Matrix> P0, int k0) = 0;


    /**
     * Calls the callable set via onResult() on all remaining smoother solutions.
     */
    virtual void endSmoother() = 0;
};


/**
 * Sequential ensemble filter interface.
 */
class ENDAS_DLL SequentialEnsembleFilter
{
public:

    /** 
     * Function called from endAnalysis() when solution is available. 
     */
    typedef std::function<void(const Ref<const Array2d> E, int k)> OnResultFn;

    SequentialEnsembleFilter();
    ~SequentialEnsembleFilter();
    SequentialEnsembleFilter(const SequentialEnsembleFilter&) = delete;


    /**
     * Sets the callable that will be called every time a filter or smoother solution is available.
     * 
     * If lag is 0, `fn` is called on the filter solutions. More usefully, when lag > 0, `fn` will be 
     * called for every smoother solution that becomes available.
     * 
     * @param fn    Any callable compatible with OnResultFn.
     * 
     * The callable is copied.
     */
    virtual void onResult(OnResultFn fn);


    /**
     * Initiates the analysis/update step.
     * 
     * The reference to the ensemble `E` is stored and will be updated as observations are being 
     * assimilated. 
     * 
     * @param E     The ensemble of state vectors at time step k.
     * @param k     The time stepping index.
     * 
     */ 
    virtual void beginAnalysis(Ref<Array2d> E, int k) = 0;

    /**
     * Assimilates observations in the vector `z` to the state vector `x`.
     * 
     * @param z     Array of observations.
     * @param H     Observation operator.
     * @param R     Observation error covariance operator.
     */
    virtual void assimilate(const Ref<const Array> z, const ObservationOperator& H, 
                            const CovarianceOperator& R) = 0;
    
    /**
     * Must be called after all observations have been assimilated.
     */ 
    virtual void endAnalysis() = 0;


protected:
    OnResultFn mOnResultFn;

};



/**
 * Sequential ensemble smoother interface.
 */
class ENDAS_DLL SequentialEnsembleSmoother : public SequentialEnsembleFilter
{
public:

    /** 
     * Must be called first to provide the initial ensemble.
     * 
     * @param E0    Initial ensemble at time step `k==0`.
     * @param k0    The initial time step index (typically 0).
     */
    virtual void beginSmoother(const Ref<const Array2d> E0, int k0) = 0;


    /**
     * Calls the callable set via onResult() on all remaining smoother solutions.
     */
    virtual void endSmoother() = 0;
};



ENDAS_DLL void ensembleForecast(Ref<Array2d> E, const GenericEvolutionModel& model,
                                const CovarianceOperator& Q, int k, double dt);




}

#endif