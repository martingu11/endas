#include <Endas/DA/Ensemble.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;



void endas::generateEnsemble(const Ref<const Array> u, const CovarianceOperator& cov, Ref<Array2d> out)
{
    ENDAS_ASSERT(out.rows() == u.size());
    cov.randomMultivariateNormal(out);
    toAnomaly(out, out); // Ensure zero mean of the sample
    out.colwise()+= u;

    /// @todo Implement oversampling and selection of orthogonal ensemble members.
}


void endas::toAnomaly(const Ref<const Array2d> E, Ref<Array2d> out)
{
    Array Eu = E.rowwise().mean();
    out = E.colwise() - Eu;
}



void endas::inflateInPlace(Ref<Array2d> E, double k)
{
    Array Eu = E.rowwise().mean();
    inflateInPlace(E, k, Eu);
}

void endas::inflateInPlace(Ref<Array2d> E, double k, Ref<Array> Eu)
{
    E = ((E.colwise() - Eu) * k).colwise() + Eu;
}

