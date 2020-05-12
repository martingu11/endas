/**
 * @file Model.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_MODEL_HPP__
#define __ENDAS_DA_MODEL_HPP__

#include <Endas/Core/Core.hpp>
#include <Endas/DA/CovarianceOperator.hpp>
#include <Endas/DA/ObservationOperator.hpp>

#include <functional>


namespace endas
{


/**
 * Base class for state evolution models.
 * 
 * The EvolutionModel and derived base classes are used primarily with the toy models included with
 * ENDAS. For use with real-world models, it is usually not necessary to wrap your model in a class 
 * derived from EvolutionModel (or in any class for that matter). 
 * 
 * For convenience, EvolutionModel is implicitly constructible from any callable (function, functor) 
 * with signature 
 * 
 *     void(Ref<Array2d> x, int k, double dt)
 * 
 * The callable object is stored via ``std::function`` and is therefore copied. If implementing your
 * own model as a class, it is preferred to derive from EvolutionModel to avoid the copy.
 */
class ENDAS_DLL EvolutionModel
{
public:

    EvolutionModel();


    /**
     * Constructs EvolutionModel as a wrapper for a callable object. 
     */
    EvolutionModel(std::function<void(Ref<Array2d> x, int k, double dt)> fn);

    virtual ~EvolutionModel();

   /**
     * Propagates state vector from time `t` to `t+dt`.
     * 
     * The state vector `x` is modified in-place. The parameter `store` informs the model whether
     * the model trajectory should be stored for linearization. Therefore, it only applies to models
     * that also implement the LinearizedEvolutionModel interface and the parameter can be ignored 
     * otherwise.
     * 
     * @param x      The state vector or ensemble of state vectors.
     * @param k      The time step index at which the model is applied.
     * @param dt     The time stepping interval.
     * @param store  Specifies whether the trajectory of the model should be stored. 
     * 
     */
    virtual void operator()(Ref<Array2d> x, int k, double dt, bool store) const = 0;

private:
    std::function<void(Ref<Array2d> x, int k, double dt)> mFn;
};



/**
 * Linearized state evolution model. 
 * 
 * The model implements both the state propagation from time step `k` to `k+1` and its linearization
 * in form of the tangent-linear and adjoint at time `k`.
 * 
 */ 
class ENDAS_DLL LinearizedEvolutionModel : public EvolutionModel
{
public:

   
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
     * Use ``std::move()`` to contructs MatrixModel without copying the passed in matrix as follows:
     * 
     *     Matrix M = ...;
     *     MatrixModel model(M);               // M is copied
     *     MatrixModel model(std::move(M));    // M is moved
     * 
     * Please note that ``M`` can not used from the calling scope after ``std::move()``. Use 
     * MatrixModel::get() to obtain reference to it.
     * 
     * @param M     Model matrix.
     */
    MatrixModel(Matrix M);

    /** Returns reference to the internal model matrix. */
    const Matrix& get() const;
    
    virtual void operator()(Ref<Array2d> x, int k, double dt, bool store = true) const override;
    virtual void tl(Ref<Array2d> x, int k) const override;
    virtual void adj(Ref<Array2d> x, int k) const override;

private:
    Matrix mModel;
};




}

#endif