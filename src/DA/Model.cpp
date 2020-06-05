#include <Endas/DA/Model.hpp>

using namespace std;
using namespace endas;


EvolutionModel::EvolutionModel()
{ }

EvolutionModel::EvolutionModel(std::function<void(Ref<Array2d> x, int k, double dt)> fn)
: mFn(fn)
{ }

EvolutionModel::~EvolutionModel()
{ }

void EvolutionModel::operator()(Ref<Array2d> x, int k, double dt, bool store) const
{
    mFn(x, k, dt);
}


MatrixModel::MatrixModel(Matrix M)
: mModel(std::move(M))
{ }

const Matrix& MatrixModel::get() const
{
    return mModel;
}

Matrix& MatrixModel::get()
{
    return mModel;
}


void MatrixModel::operator()(Ref<Array2d> x, int k, double dt, bool store) const
{
    x.matrix() = this->mModel * x.matrix();

}

void MatrixModel::tl(Ref<Array2d> x, int k) const
{
    x.matrix() = this->mModel * x.matrix();
}

void MatrixModel::adj(Ref<Array2d> x, int k) const
{
    x.matrix() = x.matrix() * this->mModel.transpose();
}

