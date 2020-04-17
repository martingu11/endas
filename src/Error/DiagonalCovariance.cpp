
#include <Endas/Error/DiagonalCovariance.hpp>
#include <Endas/Endas.hpp>
#include <Endas/Random/Random.hpp>


using namespace std;
using namespace endas;

DiagonalCovariance::DiagonalCovariance(const Ref<const Array> diag, bool isInverse)
{
    ENDAS_ASSERT(diag.size() > 0);

    mSize = diag.size();
    if (isInverse) mInvDiag = mDiag;
    else mDiag = diag;
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
    rng_engine_t& gen = getRngEngine();
    static normal_distribution<> dist;

    if (mDiagSD.size() == 0)
    {
        mDiagSD = diagonal().sqrt();
    }

    out = out.unaryExpr([&](real_t x){ return dist(gen); });
    out.colwise() *= mDiagSD;
}

void DiagonalCovariance::solve(const Ref<const Matrix> b, Ref<Matrix> out) const
{
    auto invdiag = inverseDiagonal();
    for (int i = 0; i != b.rows(); i++)
    {
        out.row(i) = invdiag(i) * b.row(i);
    }
}

void DiagonalCovariance::toMatrixView(Ref<Matrix> out) const
{
    ENDAS_ASSERT(out.rows() == mSize);
    ENDAS_ASSERT(out.cols() == mSize);

    out.fill(0);
    out.diagonal() = this->diagonal();
}


