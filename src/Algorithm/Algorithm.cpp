#include <Endas/Algorithm/Algorithm.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;


SequentialFilter::SequentialFilter()
{ }

SequentialFilter::~SequentialFilter()
{ }

void SequentialFilter::onResult(OnResultFn fn)
{
    mOnResultFn = fn;
}


SequentialEnsembleFilter::SequentialEnsembleFilter()
{ }

SequentialEnsembleFilter::~SequentialEnsembleFilter()
{ }

void SequentialEnsembleFilter::onResult(OnResultFn fn)
{
    mOnResultFn = fn;
}




void endas::ensembleForecast(Ref<Array2d> E, const GenericEvolutionModel& model,
                             const CovarianceOperator& Q, int k, double dt)
{
    ENDAS_PERF_SCOPE(EnsembleForecast);
    model(E, k, dt);

    {
        ENDAS_PERF_SCOPE(ModelEnsemblePerturbation);
        Array2d pert(E.rows(), E.cols());
        Q.randomMultivariateNormal(pert);
        E+= pert;
    }
}
