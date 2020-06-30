
#include <Endas/DA/ObservationOperator.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;


ObservationOperator::~ObservationOperator() 
{ }

bool ObservationOperator::isLinear() const { return false; }
bool ObservationOperator::isMatrix() const { return false; } 

Matrix ObservationOperator::toDenseMatrix() const
{
    ENDAS_NOT_SUPPORTED("Observation operator does not implement toDenseMatrix()");
}



shared_ptr<const ObservationOperator> ObservationOperator::subset(const IndexArray& indices) const
{
    // Fallback method using dense matrix
    if (isMatrix())
    {
        const Matrix& H = this->toDenseMatrix();

        Matrix Hsub(indices.size(), H.cols());
        selectRows(H, indices, Hsub);
        return make_shared<MatrixObservationOperator>(move(Hsub));
    }
    else
    {
        return nullptr;
    }
}



//-------------------------------------------------------------------------------------------------
// MatrixObservationOperator
//-------------------------------------------------------------------------------------------------


MatrixObservationOperator::MatrixObservationOperator(Matrix H)
: mH(move(H))
{ }

index_t MatrixObservationOperator::nobs() const { return mH.rows(); }
index_t MatrixObservationOperator::nstate() const { return mH.cols(); }

bool MatrixObservationOperator::isLinear() const { return true; }
bool MatrixObservationOperator::isMatrix() const { return true; }

void MatrixObservationOperator::apply(const Ref<const Array2d> x, Ref<Array2d> out) const
{
    out.matrix().noalias() = this->mH * x.matrix();
}

Matrix MatrixObservationOperator::toDenseMatrix() const
{
    return this->mH;
}


//-------------------------------------------------------------------------------------------------
// CustomObservationOperator
//-------------------------------------------------------------------------------------------------



CustomObservationOperator::CustomObservationOperator(int nobs, int nstate, bool isLinear, 
                                                     const ApplyFn& fn)
: mNObs(nobs), mNState(nstate), mIsLinear(isLinear), mApplyFn(fn)
{ }


index_t CustomObservationOperator::nobs() const { return mNObs; }
index_t CustomObservationOperator::nstate() const { return mNState; }

bool CustomObservationOperator::isLinear() const { return mIsLinear; }
bool CustomObservationOperator::isMatrix() const { return false; }

void CustomObservationOperator::apply(const Ref<const Array2d> x, Ref<Array2d> out) const
{
    mApplyFn(x, out);
}




