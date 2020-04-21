
#include <Endas/DA/ObservationOperator.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;

ObservationOperator::~ObservationOperator() 
{ }




MatrixObservationOperator::MatrixObservationOperator(const Ref<const Matrix> H)
: mH(H)
{ }

MatrixObservationOperator::MatrixObservationOperator(int rows, int cols, 
                                                     const function<void(Ref<Matrix>)>& H)
: mH(rows, cols)
{
    H(this->mH);
}

int MatrixObservationOperator::nobs() const { return mH.rows(); }

bool MatrixObservationOperator::isLinear() const
{
    return true;
}

bool MatrixObservationOperator::isMatrix() const
{
    return true;    
}

void MatrixObservationOperator::apply(const Ref<const Array2d> x, int k, Ref<Array2d> out) const
{
    out.matrix().noalias() = this->mH * x.matrix();
}

void MatrixObservationOperator::toMatrix(Matrix& out) const
{
    out = this->mH;
}

