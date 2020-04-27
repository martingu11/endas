
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
    RandomNumberGenerator& gen = getRandomNumberGenerator();

    for (int i = 0; i != out.cols(); i++)
    {
        auto col = out.col(i);

        gen.standardNormal(col);
        col.matrix().noalias() = mTransform * mX.matrix();
        if (mMean.size() > 0) col+= mMean;
    }
}

