/**
 * @file DiagonalCovariance.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_DIAGONAL_COVARIANCE_OPERATOR_HPP__
#define __ENDAS_DIAGONAL_COVARIANCE_OPERATOR_HPP__

#include "CovarianceOperator.hpp"

namespace endas
{


/**
 * Implements diagonal covariance matrix.
 *
 * The covariance operator is internally represented by the diagonal array. Currently only the main 
 * diagonal is supported but other diagonals could be added, if needed. The operator supports all
 * methods of CovarianceOperator. The covariance can be instantiated with either the array of diagonal
 * elements or the reciprocal (inverse) array. This can prevent numerical issues in situations where 
 * the inverse coefficients are near zero (thus leading to very large coefficients on the original 
 * diagonal) and if only the inverse coefficients are needed, such as when only 
 * CovarianceOperator::solve() is called.
*/
class ENDAS_DLL DiagonalCovariance : public CovarianceOperator
{
public:

    /**
     * Constructs diagonal covariance with given diagonal or its inverse.
     * The passed array is copied.
     * 
     * @param diag       Array of diagonal elements.
     * @param isInverse  If `true`, the `diag` array is assumed to hold the inverse diagonal 
     *                   coefficients, otherwise the actual coefficients are assumed.
     */
    DiagonalCovariance(const Ref<const Array> diag, bool isInverse = false);

    /**
     * Returns reference to the matrix diagonal.
     *
     * The diagonal is constructed on demand (if DiagonalCovariance was initialized with 
     * inverse diagonal).
     */ 
    const Ref<const Array> diagonal() const;

    /**
     * Returns reference to the matrix diagonal.
     * 
     * The inverse diagonal is constructed on demand (if DiagonalCovariance was initialized with 
     * actual diagonal).
     */ 
    const Ref<const Array> inverseDiagonal() const;

    virtual int size() const override;
    virtual bool isDiagonal() const override;
    virtual bool mcOnly() const override;
    virtual void randomMultivariateNormal(Ref<Array2d> out) const override;
    virtual void solve(const Ref<const Matrix> b, Ref<Matrix> out) const override;
    virtual void toMatrixView(Ref<Matrix> out) const override;

private:
   
    int mSize;
    mutable Array mDiag;
    mutable Array mInvDiag;
    mutable Array mDiagSD;
};




}

#endif