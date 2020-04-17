
#include <Endas/Error/CovarianceOperator.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;

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

void CovarianceOperator::toMatrix(Matrix& out) const
{
    out.resize(this->size(), this->size());
    this->toMatrixView(out);
}

void CovarianceOperator::toMatrixView(Ref<Matrix> out) const
{
    ENDAS_NOT_SUPPORTED("Covariance operator does not support toMatrix()");
}

void CovarianceOperator::solve(const Ref<const Matrix> b, Ref<Matrix> out) const
{
    ENDAS_NOT_SUPPORTED("Covariance operator does not support solve()");
}

