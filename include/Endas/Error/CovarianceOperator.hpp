/**
 * @file CovarianceOperator.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_COVARIANCE_OPERATOR_HPP__
#define __ENDAS_COVARIANCE_OPERATOR_HPP__

#include <Endas/Core/LinAlg.hpp>


namespace endas
{

/**
 * Abstract representation of a covariance matrix.
 *
 * The base class defines interface required by covariance matrix implementations. 
 */
class ENDAS_DLL CovarianceOperator
{
public:

    virtual ~CovarianceOperator();

    /**
     * Returns the size of the space the covariance represents. 
     */
    virtual int size() const = 0;

    /**
     * Returns `true` if the covariance operator can be represented by a diagonal matrix.
     * The default implementation returns `false`.
     */
    virtual bool isDiagonal() const;

    /**
     * Returns `true` if this covariance operator only supports Monte-Carlo sampling via the 
     * randomMultivariateNormal() call. 
     */
    virtual bool mcOnly() const;

    /**
     * Generates random sample from a multivariate Normal distribution with zero mean and covariance
     * given by this CovarianceMatrix instance.
     * 
     * @param N     The number of independent samples to draw.
     * @param out   Pre-allocated array instance where the sample is stored.
     * 
     * The `out` matrix must be of size (n, N), where `n = this->size()`.
     */
    virtual void randomMultivariateNormal(Ref<Array2d> out) const = 0;


    virtual void solve(const Ref<const Matrix> b, Ref<Matrix> out) const;


    /**
     * Computes `x = x + R * mult`, where `R` is the covariance matrix represented by this operator.
     * 
     * Please note that not all covariance operator implementations support this. Use `mcOnly()` to 
     * check if the matrix form is available. 
     * 
     * @param x   Reference to the array to modify.
     * @throw NotSupportedError if the operation is not supported.
     */
    virtual void addTo(Ref<Array2d> x, double mult = 1.0) const;


    /**
     * Returns dense matrix representation of the covariance. 
     * 
     * Please note that not all covariance operator implementations support this. Use `mcOnly()` to 
     * check if the matrix form is available. 
     * 
     * @param out   Reference to a matrix expression where the data is to be written  
     * @throw NotSupportedError if the operation is not supported.
     * 
     * @note: This method can write the covariance matrix to any matrix or array expression that is 
     * of the correct size (size() x size()). Passing wrongly-sized matrix expression as `out` will 
     * trigger runtime assertion. Use toMatrix() to have the matrix resized automatically.
     */
    virtual void toMatrixView(Ref<Matrix> out) const;

    /**
     * Returns dense matrix representation of the covariance. 
     * 
     * Unlike toMatrixView(), toMatrix() can only write to actual Matrix instances, not expressions.
     * The matrix is resized automatically to the correct size.
     */ 
    void toMatrix(Matrix& out) const;


};




}

#endif