/**
 * @file ObservationOperator.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_DA_OBSERVATION_OPERATOR_HPP__
#define __ENDAS_DA_OBSERVATION_OPERATOR_HPP__

#include <Endas/Core/LinAlg.hpp>


namespace endas
{

/** 
 * @addtogroup da
 * @{ 
 */


/**
 * Abstract base class for observation operators.
 *
 * Observation operators provide mapping from the state space to the observation space. In other
 * words, the operator transforms the state vector to the corresponding set of observations that
 * we would expect.
 */
class ENDAS_DLL ObservationOperator
{
public:

    virtual ~ObservationOperator();

    /** Returns the size of the observation space spanned by this operator. */
    virtual index_t nobs() const = 0;

    /** Returns the size of the state space spanned by this operator. */
    virtual index_t nstate() const = 0;

    /** 
     * Returns `true` if the operator is linear. 
     * Default implementation returns `false`. 
     */
    virtual bool isLinear() const; 

    /** 
     * Returns `true` if the operator can be represented by a matrix. 
     * Default implementation returns `false`.
     */
    virtual bool isMatrix() const;   

    /**
     * Applies the observation operator to array or matrix `x`.
     * The array or matrix `x` is updated in-place, i.e. this computes `x = H(x)`. 
     * 
     * @param x     The array or matrix to apply the operator to.
     * @param k     Current time step index.
     * @param out   Pre-alloccated array of size nobs() x `x.cols()` where the result is to be stored.
     */
    virtual void apply(const Ref<const Array2d> x, int k, Ref<Array2d> out) const = 0;


    /**
     * Returns dense matrix representation of the operator. 
     * 
     * Please note that not all observation operator implementations support this. Use isMatrix() 
     * to check if asMatrix() is allowed.
     * 
     * @return  Reference to a matrix
     * @throw   NotSupportedError if the operation is not supported.
     */
    virtual const Matrix& asMatrix() const;


    /**
     * Returns observation operator that is a subset of this operator.
     * 
     * Please note that not all observation operator implementations support this. If isMatrix() 
     * returns `true`, subsetting will be supported. Otherwise subset() may return `nullptr`.
     * 
     * @param indices  Array of indices that select which observations are included in the subset.
     * 
     * @note Depending on the implementation, the returned instance _may_ share data with this 
     * operator. Therefore, it is advised to ensure both operators remain alive long enough in the
     * calling code.
     *  
     * @see The ObservationManager interface offers more flexibility for handling localization of
     * observations and should generally be preferred for more complex cases.
     */
    virtual std::shared_ptr<const ObservationOperator> subset(const IndexArray& indices) const;

};



/**
 * Simple observation operator represented by a matrix.
 */
class ENDAS_DLL MatrixObservationOperator : public ObservationOperator
{
public:

    /** 
     * MatrixObservationOperator constructor.
     * 
     * @param H     Operator matrix. The referenced matrix instance is copied or moved.
     */
    MatrixObservationOperator(Matrix H);


    virtual index_t nobs() const override;
    virtual index_t nstate() const override;

    virtual bool isLinear() const override; 
    virtual bool isMatrix() const override;   

    virtual void apply(const Ref<const Array2d> x, int k, Ref<Array2d> out) const override;
    virtual const Matrix& asMatrix() const override;

private:
    Matrix mH;
};




/**
 * Observation operator represented by a user-defined callable.
 */
class ENDAS_DLL CustomObservationOperator : public ObservationOperator
{
public:

    typedef std::function<void(const Ref<const Array2d> x, int k, Ref<Array2d> out)> ApplyFn;

    /** 
     * MatrixObservationOperator constructor.
     * 
     * @param H     Operator matrix. The referenced matrix instance is copied.
     */
    CustomObservationOperator(int nobs, int nstate, bool isLinear, const ApplyFn& fn);


    virtual index_t nobs() const override;
    virtual index_t nstate() const override;

    virtual bool isLinear() const override; 
    virtual bool isMatrix() const override;   

    virtual void apply(const Ref<const Array2d> x, int k, Ref<Array2d> out) const override;

private:
    index_t mNObs, mNState;
    bool mIsLinear;
    ApplyFn mApplyFn;
};




/** @} */

}

#endif