
#include <Endas/Random/Random.hpp>
#include <Endas/Random/MultivariateRandomNormal.hpp>

using namespace std;
using namespace endas;


MultivariateRandomNormal::MultivariateRandomNormal(const Ref<const Matrix> cov)
: mX(cov.rows())
{ 
    Eigen::SelfAdjointEigenSolver<Matrix> eigenSolver(cov);
    mTransform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
}

MultivariateRandomNormal::MultivariateRandomNormal(const Ref<const Array> mean, 
                                                   const Ref<const Matrix> cov)
: MultivariateRandomNormal(cov)
{ 
    mMean = mean;
}

void MultivariateRandomNormal::operator()(Ref<Array2d> out) const
{
    rng_engine_t& gen = getRngEngine();
    static normal_distribution<> dist;

    for (int i = 0; i != out.cols(); i++)
    {
        auto col = out.col(i);

        mX = mX.unaryExpr([&](real_t x) { return dist(gen); });
        col.matrix().noalias() = mTransform * mX.matrix();
        if (mMean.size() > 0) col+= mMean;
    }
}

