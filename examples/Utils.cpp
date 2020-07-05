#include "Utils.hpp"
#include <Endas/Core/Ensemble.hpp>
#include <Endas/DA/ObservationOperator.hpp>
#include <Endas/DA/CovarianceOperator.hpp>
#include <Endas/DA/SimpleObservationManager.hpp>

using namespace std;
using namespace endas;


tuple<Array2d, vector<Array>> 
endas::generateExampleData(int numSteps, const Ref<const Array> x0, 
                           const EvolutionModel& model, double dt,
                           const CovarianceOperator& Q, 
                           const ObservationOperator& H, 
                           const CovarianceOperator& R,
                           int numSpinupSteps, int obsInterval)
{
    ObservationOpAndCov HR = 
    {
        shared_ptr_wrap(H), shared_ptr_wrap(R), 
    };

    return generateExampleData(
        numSteps, x0, model, dt, Q,
        [&](int k) { return HR; },
        numSpinupSteps, obsInterval
    );
}



std::tuple<Array2d, std::vector<Array>> 
endas::generateExampleData(int numSteps, const Ref<const Array> x0,
                           const EvolutionModel& model, double dt,
                           const CovarianceOperator& Q, 
                           std::function<ObservationOpAndCov(int k)> HRfun,
                           int numSpinupSteps, int obsInterval)
{
    int n = x0.size();
    //int nobsSteps = (int)ceil((numSteps-1) / obsInterval);

    ENDAS_ASSERT(n > 0);
    ENDAS_ASSERT(numSteps > 0);

    // The 'true' system state evolution
    Array2d xtAll = Array2d::Zero(n, numSteps);
    xtAll.col(0) = x0;

    // The observation vectors
    vector<Array> zAll;

    Array x = x0;
    Array xnoise(n);
    
    // Model spin-up
    for (int k = 0; k != numSpinupSteps; k++)
    {
        model(x, k, dt, false);
    }

    // At k==0 we have no observations
    zAll.push_back(Array());

    // Generate true state and observations
    for (int k = 1; k != numSteps; k++)
    {
        model(x, k, dt, false);

        Q.randomMultivariateNormal(xnoise);
        x+= xnoise;
        xtAll.col(k) = x;

        // Note: Generate noise even if we don't use it for deterministic results regardless of 
        // obsInterval

        auto HR = HRfun(k);
        ENDAS_ASSERT(HR.H);
        ENDAS_ASSERT(HR.R);

        Array znoise(HR.H->nobs());
        HR.R->randomMultivariateNormal(znoise);

        if (k % obsInterval == 0)
        {
            Array z(HR.H->nobs());
            HR.H->apply(x, z);
            z += znoise;
            zAll.push_back(move(z));
        }
        else 
            zAll.push_back(Array());
    }

    //ENDAS_ASSERT(obsSteps.size() == zAll.cols());
    return make_tuple(xtAll, zAll);
}

