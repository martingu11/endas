/**
 * @file CovarianceOperator.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_COVARIANCE_OPERATOR_HPP__
#define __ENDAS_DA_COVARIANCE_OPERATOR_HPP__

#include <Endas/Core/LinAlg.hpp>
#include <memory>


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
     * Computes fused multiply-add (`A = A + B * c)`, where `B` is the covariance matrix 
     * represented by this operator.
     * 
     * Please note that not all covariance operator implementations support this. Use `mcOnly()` to 
     * check if the matrix form is available. 
     * 
     * @param A     Reference to the array to add to.
     * @param c     Coefficient-wise multiplication factor applied to `B`.
     * 
     * @throw NotSupportedError if the operation is not supported.
     */
    virtual void fmadd(Ref<Array2d> A, double c = 1.0) const;


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
    //virtual void toMatrixView(Ref<Matrix> out) const;

    /**
     * Returns dense matrix representation of the covariance operator. 
     */ 
    virtual const Matrix& asMatrix() const;


    /**
     * Returns covariance operator that represents covariance of a subset of the original space.
     * 
     * Please note that not all observation operator implementations support this. If mcOnly() 
     * returns `false`, subsetting will be supported. Otherwise subset() may return `nullptr`.
     * 
     * @param indices  Array of indices that select which elements are included in the subset.
     * 
     * @note Depending on the implementation, the returned instance _may_ share data with this 
     * operator. Therefore, it is advised to ensure both operators remain alive long enough in the
     * calling code.
     *  
     * @see The ObservationManager interface offers more flexibility for handling localization of
     * covariance and should generally be preferred for more complex cases.
     */
    virtual std::shared_ptr<const CovarianceOperator> subset(const Ref<const Array> indices) const;

};


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
    virtual void fmadd(Ref<Array2d> A, double c = 1.0) const override;
    virtual const Matrix& asMatrix() const override;

    virtual std::shared_ptr<const CovarianceOperator> subset(const Ref<const Array> indices) const override;

private:
   
    bool mInitWithInverse;
    int mSize;
    mutable Array mDiag;
    mutable Array mInvDiag;
    mutable Array mDiagSD;
    mutable Matrix mDense; // Constructed when needed
};



/**
 * Implements dense covariance matrix.
 *
 * Use only on small covariances!
 */
class ENDAS_DLL DenseCovariance : public CovarianceOperator
{
public:

    /**
     * Constructs diagonal covariance with given diagonal or its inverse.
     * The passed array is copied.
     * 
     * @param P       Covariance matrix. The matrix is copied.
     */
    DenseCovariance(const Ref<const Matrix>& P);

    DenseCovariance(const Matrix&& P);

    ~DenseCovariance();

    virtual int size() const override;
    virtual bool isDiagonal() const override;
    virtual bool mcOnly() const override;
    virtual void randomMultivariateNormal(Ref<Array2d> out) const override;
    virtual void solve(const Ref<const Matrix> b, Ref<Matrix> out) const override;

    virtual void fmadd(Ref<Array2d> A, double c = 1.0) const override;
    virtual const Matrix& asMatrix() const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};







}

#endif