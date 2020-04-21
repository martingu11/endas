#include <Endas/DA/Model.hpp>

using namespace std;
using namespace endas;



MatrixModel::MatrixModel(const Ref<const Matrix> M)
: mModel(M)
{ }

MatrixModel::MatrixModel(int n, const std::function<void(Ref<Matrix>)> M)
: mModel(n, n)
{
    M(this->mModel);
}

const Matrix& MatrixModel::get() const
{
    return mModel;
}

void MatrixModel::apply(Ref<Array2d> x, int k, double dt, bool store) const
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

