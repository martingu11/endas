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
    virtual int nobs() const = 0;

    /** Returns `true` if the operator is linear. */
    virtual bool isLinear() const = 0; 

    /** Returns `true` if the operator can be represented by a matrix. */
    virtual bool isMatrix() const = 0;   

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
     * to check if toMatrix() is allowed.
     * 
     * @param out   Reference to a matrix expression where the data is to be written  
     * @throw NotSupportedError if the operation is not supported.
     */
    virtual void toMatrix(Matrix& out) const = 0;

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
     * @param H     Operator matrix. The referenced matrix instance is copied.
     */
    MatrixObservationOperator(const Ref<const Matrix> H);


    /**
     * MatrixObservationOperator constructor.
     * 
     * @param rows  Number of rows.
     * @param cols  Number of columns.
     * @param H     Any callable with signature `void(Ref<Matrix>)` that populates the operator
     *              coefficients.
     */
    MatrixObservationOperator(int rows, int cols, const std::function<void(Ref<Matrix>)>& M);


    virtual int nobs() const override;
    virtual bool isLinear() const override; 
    virtual bool isMatrix() const override;   

    virtual void apply(const Ref<const Array2d> x, int k, Ref<Array2d> out) const override;
    virtual void toMatrix(Matrix& out) const override;

private:

    Matrix mH;
};




}

#endif