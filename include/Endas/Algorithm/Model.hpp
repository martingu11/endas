/**
 * @file Model.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_ALGORITHM_MODEL_HPP__
#define __ENDAS_ALGORITHM_MODEL_HPP__

#include <Endas/Core/Core.hpp>
#include <Endas/Core/LinAlg.hpp>
#include <Endas/Error/CovarianceOperator.hpp>
#include <Endas/Observation/ObservationOperator.hpp>

#include <functional>


namespace endas
{


/**
 * Generic state evolution model.
 * 
 * The model can be any callable with signature `void(x, k, dt)` that updates the state vector
 * or ensemble `x` from time step `k` to `k+1`. 
 * 
 * @param x     Column vector or matrix holding the state vector or ensemble of state vectors, 
 *              respectively.
 * @param k     Time step index corresponding to `x` before the model call.
 * @param dt    Time stepping interval to be applied.
 */ 
typedef std::function<void(Ref<Array2d> x, int k, double dt)> GenericEvolutionModel; 


/**
 * Linearized state evolution model. 
 * 
 * The model implements both the state propagation from time step `k` to `k+1` and its linearization
 * in form of the tangent-linear and adjoint at time `k`.
 * 
 */ 
class LinearizedEvolutionModel
{
public:

    virtual ~LinearizedEvolutionModel() { }

    /**
     * Propagates state vector from time `t` to `t+dt`.
     * The state vector `x` is modified in-place.
     * 
     * @param x     The state vector or ensemble of state vectors.
     * @param k     The time step index at which the model is applied.
     * @param dt    The time stepping interval.
     * @param store Specifies whether data needed for tl() and adj() should be stored.
     *              Can be set to `false` if neither tl() nor adj() is used.
     */
    virtual void apply(Ref<Array2d> x, int k, double dt, bool store = true) const = 0;


    /**
     * Applies tangent-linear of the model at step `k` to `x`.
     * 
     * @param x     Matrix to apply the tangent-linear to.
     * @param k     The time step index at which the tangent-linear is applied.
     */
    virtual void tl(Ref<Array2d> x, int k) const = 0;

    
    /**
     * Applies adjoint of the model at step `k` to `x`.
     * 
     * @param x     Matrix to apply the adjoint to.
     * @param k     The time step index at which the adjoint is applied.
     */
    virtual void adj(Ref<Array2d> x, int k) const = 0;


    /**
     * Notifies the model that trajectory data for time step `k` is no longer needed.
     * 
     * @param k     The time step index for which data can be dropped.
     * 
     * The default implementation does nothing.
     */ 
    virtual void stepFinished(int k) const
    { }

    // Implicit conversion to GenericEvolutionModel. 
    operator GenericEvolutionModel() const
    {
        return std::bind(&LinearizedEvolutionModel::apply, this, 
                         std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 
                         false);
    }


};



/**
 * Trivial state evolution model represented by a matrix.
 */ 
class ENDAS_DLL MatrixModel : public LinearizedEvolutionModel
{
public:

    /** 
     * MatrixModel constructor.
     * 
     * @param M     Model matrix. The referenced matrix instance is copied.
     */
    MatrixModel(const Ref<const Matrix> M);

    /**
     * MatrixModel constructor.
     * 
     * @param n     State size.
     * @param M     Any callable with signature `void(Ref<Matrix>)` that populates the model
     *              matrix coefficients.
     */
    MatrixModel(int n, const std::function<void(Ref<Matrix>)> M);

    /** Returns reference to the internal model matrix. */
    const Matrix& get() const;
    
    virtual void apply(Ref<Array2d> x, int k, double dt, bool store = true) const override;
    virtual void tl(Ref<Array2d> x, int k) const override;
    virtual void adj(Ref<Array2d> x, int k) const override;

private:
    Matrix mModel;
};




}

#endif