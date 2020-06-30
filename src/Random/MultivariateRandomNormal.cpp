
#include <Endas/Random/Random.hpp>
#include <Endas/Random/MultivariateRandomNormal.hpp>
#include "../Compatibility.hpp"

#include <iostream>

using namespace std;
using namespace endas;


struct MultivariateRandomNormal::Data
{
    Matrix cov;
    Array mean;

    // Declare after `cov` as it uses Ref to it!
    Eigen::LLT<Ref<Matrix>> llt;

    Data(Matrix _cov)
    : cov(move(_cov)),
      llt(this->cov)
    { }
};


MultivariateRandomNormal::MultivariateRandomNormal(Matrix cov)
: mData(make_unique<Data>(move(cov)))
{ 
    //mData->cov = move(cov);
    //mData->llt.compute(mData->cov);
    //mL = llt.matrixL();
    //Eigen::SelfAdjointEigenSolver<Matrix> eigenSolver(cov);
    //mTransform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
}

MultivariateRandomNormal::MultivariateRandomNormal(Array mean, Matrix cov)
: MultivariateRandomNormal(move(cov))
{ 
    mData->mean = move(mean);
}

MultivariateRandomNormal::~MultivariateRandomNormal()
{ }


void MultivariateRandomNormal::sample(Ref<Array2d> out) const
{
    RandomNumberGenerator& gen = getRandomNumberGenerator();

    Matrix X(mData->cov.rows(), 1);
    for (int i = 0; i != out.cols(); i++)
    {
        gen.standardNormal(X);
        out.col(i).matrix().noalias() = mData->llt.matrixL() * X;
        //col.matrix().noalias() = mTransform * mX.matrix();
        if (mData->mean.size() > 0) out.col(i)+= mData->mean;
    }
}

