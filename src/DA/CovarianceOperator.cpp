
#include <Endas/DA/CovarianceOperator.hpp>
#include <Endas/Endas.hpp>
#include <Endas/Random/Random.hpp>
#include <Endas/Random/MultivariateRandomNormal.hpp>

#include "../Compatibility.hpp"


using namespace std;
using namespace endas;


//-----------------------------------------------------------------------------
// DenseCovariance
//-----------------------------------------------------------------------------


CovarianceOperator::~CovarianceOperator() 
{ }

bool CovarianceOperator::isDiagonal() const
{
    return false;
}

bool CovarianceOperator::mcOnly() const
{
    return true;    
}


const Matrix& CovarianceOperator::asMatrix() const
{
    ENDAS_NOT_SUPPORTED("Covariance operator does not support asMatrix()");
}


void CovarianceOperator::solve(const Ref<const Matrix> b, Ref<Matrix> out) const
{
    ENDAS_NOT_SUPPORTED("Covariance operator does not support solve()");
}


void CovarianceOperator::fmadd(Ref<Array2d> A, double c) const
{
    ENDAS_NOT_SUPPORTED("Covariance operator does not support fmadd()");
}


shared_ptr<const CovarianceOperator> CovarianceOperator::subset(const IndexArray& indices) const
{
    // Fallback method using dense matrix
    if (!mcOnly())
    {
        const Matrix& P = this->asMatrix();
        
        Matrix Psub(indices.size(), indices.size());
        selectRowsCols(P, indices, indices, Psub);
        return make_shared<DenseCovariance>(move(Psub));
    }
    else
    {
        return nullptr;
    }
}


//-----------------------------------------------------------------------------
// DiagonalCovariance
//-----------------------------------------------------------------------------

DiagonalCovariance::DiagonalCovariance(const Ref<const Array> diag, bool isInverse)
{
    ENDAS_ASSERT(diag.size() > 0);

    mSize = diag.size();
    if (isInverse) mInvDiag = mDiag;
    else mDiag = diag;

    mInitWithInverse = isInverse;
}

int DiagonalCovariance::size() const 
{
    return mSize;
}

bool DiagonalCovariance::isDiagonal() const
{
    return true;
}

bool DiagonalCovariance::mcOnly() const
{
    return false;
}

const Ref<const Array> DiagonalCovariance::diagonal() const
{
    if (mDiag.size() == 0)
    {
        mDiag = mInvDiag.inverse();
    }
    return mDiag;
}

const Ref<const Array> DiagonalCovariance::inverseDiagonal() const
{
    if (mInvDiag.size() == 0)
    {
        mInvDiag = mDiag.inverse();
    }
    return mInvDiag;
}


void DiagonalCovariance::randomMultivariateNormal(Ref<Array2d> out) const
{
    RandomNumberGenerator& gen = getRandomNumberGenerator();

    if (mDiagSD.size() == 0)
    {
        mDiagSD = diagonal().sqrt();
    }
    gen.standardNormal(out);
    out.colwise() *= mDiagSD;
}

void DiagonalCovariance::solve(const Ref<const Matrix> b, Ref<Matrix> out) const
{
    auto invdiag = inverseDiagonal();
    diagXdense(invdiag, b, out);
}


void DiagonalCovariance::fmadd(Ref<Array2d> A, double c) const
{
    ENDAS_ASSERT(A.rows() == this->mSize && A.cols() == this->mSize);

    if (c == 1.0)
    {
        A.matrix().diagonal().array() += this->diagonal();
    }
    else
    {
        A.matrix().diagonal().array() += this->diagonal() * c;
    }
}

const Matrix& DiagonalCovariance::asMatrix() const
{
    if (mDense.size() == 0)
    {
        mDense = Matrix::Zero(this->size(), this->size());
        mDense.diagonal() = this->diagonal();
    }
    return mDense;
}

shared_ptr<const CovarianceOperator> DiagonalCovariance::subset(const IndexArray& indices) const
{
    Array diagSub(indices.size());
    select((mInitWithInverse)? mInvDiag : mDiag, indices, diagSub);
    
    return make_shared<DiagonalCovariance>(diagSub, mInitWithInverse);
}


//-----------------------------------------------------------------------------
// DenseCovariance
//-----------------------------------------------------------------------------

struct DenseCovariance::Data
{
    Matrix P; 

    bool haveCholP;
    Eigen::LLT<Matrix> cholP;
    std::unique_ptr<MultivariateRandomNormal> mrn;

    Data(Matrix _P) : P(move(_P)), haveCholP(false) { }
};


DenseCovariance::DenseCovariance(Matrix P)
: mData(make_unique<Data>(move(P)))
{ 
    ENDAS_ASSERT(mData->P.cols() == mData->P.rows());
}


DenseCovariance::~DenseCovariance()
{ }

int DenseCovariance::size() const { return mData->P.cols(); }
bool DenseCovariance::isDiagonal() const { return false; }
bool DenseCovariance::mcOnly() const { return false; }

void DenseCovariance::randomMultivariateNormal(Ref<Array2d> out) const
{
    if (!mData->mrn)
    {
        mData->mrn = make_unique<MultivariateRandomNormal>(mData->P);
    }
    mData->mrn->operator()(out);
}

void DenseCovariance::solve(const Ref<const Matrix> b, Ref<Matrix> out) const
{
    if (!mData->haveCholP)
    {
        /// @todo  We could do decomposition in place and reconstruct matrix in asMatrix()?
        mData->cholP.compute(mData->P);
        mData->haveCholP = true;
    }

    /// @todo solveInPlace() if b==out?
    out = mData->cholP.solve(b);
}

void DenseCovariance::fmadd(Ref<Array2d> A, double c) const
{
    ENDAS_ASSERT(A.rows() == mData->P.rows() && A.cols() == mData->P.cols());

    if (c == 1.0)
    {
        A += mData->P.array();
    }
    else
    {
        A += mData->P.array() * c;
    }
}

const Matrix& DenseCovariance::asMatrix() const
{
    return mData->P;
}



