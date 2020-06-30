
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


Matrix CovarianceOperator::toDenseMatrix() const
{
    ENDAS_NOT_SUPPORTED("Covariance operator does not support toDenseMatrix()");
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
        const Matrix& P = this->toDenseMatrix();
        
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
    out.noalias() = invdiag.matrix().asDiagonal() * b;
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

Matrix DiagonalCovariance::toDenseMatrix() const
{
    Matrix dense = Matrix::Zero(this->size(), this->size());
    dense.diagonal() = this->diagonal();
    return dense;
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

    bool haveLLT;
    Eigen::LLT<Matrix> LLT;

    Data(Matrix _P) : P(move(_P)), haveLLT(false) { }
};


DenseCovariance::DenseCovariance(Matrix P)
: mData(make_unique<Data>(move(P)))
{ 
    ENDAS_ASSERT(mData->P.cols() == mData->P.rows());
}


DenseCovariance::DenseCovariance(const SpatialStateSpace& space, const CovarianceFn& covFn,
                                 double epsilon)
: mData(make_unique<Data>(Matrix::Zero(space.size(), space.size())))
{
    Matrix& P = mData->P;
    int dim = space.dim();
    index_t n = space.size();
    const CoordinateSystem& crs = space.crs();

    // Get the coordinates of the state vector elements
    Array2d coords(dim, n);
    space.getCoords(coords);

    // Evaluate covFn(). We only need to initialize the lower triangular part of the matrix,
    // upper triangular is not used

    const IsotropicCovarianceFn* isoCovFn = dynamic_cast<const IsotropicCovarianceFn*>(&covFn);

    // If covFn() is IsotropicCovarianceFn
    if (isoCovFn)
    {
        Array dist(n);
        for (index_t i = 0; i != n; i++)
        {
            auto d = dist.head(n-i);
            crs.distance(coords.block(0, i, dim, n-i), coords.col(i), dist.head(n-i));
            isoCovFn->values(d, d);
            P.block(i, i, n-i, 1) = d;
        }
    }
    else
    {
        ENDAS_NOT_IMPLEMENTED;
    }

    // Add a small amount to the diagonal to increase chances the covariance is positive
    // definite
    if (epsilon > 0)
    {
        P.diagonal() = P.diagonal().array() += epsilon;
    }

}


DenseCovariance::~DenseCovariance()
{ }

int DenseCovariance::size() const { return mData->P.cols(); }
bool DenseCovariance::isDiagonal() const { return false; }
bool DenseCovariance::mcOnly() const { return false; }

void DenseCovariance::randomMultivariateNormal(Ref<Array2d> out) const
{
    int n = mData->P.rows();
    ENDAS_ASSERT(out.rows() == n);

    if (!mData->haveLLT)
    {
        mData->LLT.compute(mData->P);
        mData->haveLLT = true;
    }

    RandomNumberGenerator& gen = getRandomNumberGenerator();
    Matrix X(n, 1);
    
    for (int i = 0; i != out.cols(); i++)
    {
        gen.standardNormal(X);
        out.col(i).matrix().noalias() = mData->LLT.matrixL() * X;
    }
}

void DenseCovariance::solve(const Ref<const Matrix> b, Ref<Matrix> out) const
{
    if (!mData->haveLLT)
    {
        mData->LLT.compute(mData->P);
        mData->haveLLT = true;
    }

    /// @todo solveInPlace() if b==out?
    out = mData->LLT.solve(b);
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

Matrix DenseCovariance::toDenseMatrix() const
{
    // P might only have valid values in the lower triangular. Reconstruct full matrix 
    // from that.
    Matrix dense = mData->P.triangularView<Eigen::Lower>().selfadjointView();
    return move(dense);

    //return mData->P;
}




//-----------------------------------------------------------------------------
// ZeroCovariance
//-----------------------------------------------------------------------------

ZeroCovariance::ZeroCovariance(index_t size)
: mSize(size)
{ }


int ZeroCovariance::size() const { return mSize; }
bool ZeroCovariance::isDiagonal() const { return true; }
bool ZeroCovariance::mcOnly() const { return false; }

void ZeroCovariance::randomMultivariateNormal(Ref<Array2d> out) const
{
    out.fill(0);
}

void ZeroCovariance::solve(const Ref<const Matrix> b, Ref<Matrix> out) const
{
    ENDAS_NOT_SUPPORTED("solve() not supported for ZeroCovariance");
}

void ZeroCovariance::fmadd(Ref<Array2d> A, double c) const
{
    // Nothing to do
}

Matrix ZeroCovariance::toDenseMatrix() const
{
    return Matrix::Zero(mSize, mSize);
}





