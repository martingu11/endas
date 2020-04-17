#include "Utils.hpp"

using namespace std;
using namespace endas;


tuple<Array2d, Array2d, vector<int>> 
endas::generateTestData(int nsteps, const Ref<const Array> x0, 
                        const GenericEvolutionModel& model, double dt,
                        const Ref<const Matrix> H, const CovarianceOperator& Q, 
                        const CovarianceOperator& R,
                        int obsInterval)
{
    int n = x0.size();
    int nobs = H.rows();
    int nobsSteps = (int)ceil((nsteps-1) / obsInterval);

    ENDAS_ASSERT(n > 0);
    ENDAS_ASSERT(nsteps > 0);

    // The 'true' system state evolution
    Array2d xtAll = Array2d::Zero(n, nsteps);
    xtAll.col(0) = x0;

    // Array of observations (we have nobs observations at each step 1..nsteps)
    Array2d zAll = Array2d::Zero(nobs, nobsSteps);  

    // The observation times (steps)
    vector<int> obsSteps;

    Array x = x0;
    Array xnoise(n);
    Array znoise(nobs);
    for (int k = 1; k != nsteps; k++)
    {
        model(x, k, dt);
        
        Q.randomMultivariateNormal(xnoise);
        x+= xnoise;
        xtAll.col(k) = x;

        // Note: Generate noise even if we don't use it for deterministic results regardless of 
        // obsInterval
        R.randomMultivariateNormal(znoise);

        if (k % obsInterval == 0)
        {
            Array z = H * x.matrix();
            
            z += znoise;
            zAll.col(obsSteps.size()) = z;
            obsSteps.push_back(k);
        }
    }

    ENDAS_ASSERT(obsSteps.size() == zAll.cols());
    return make_tuple(xtAll, zAll, obsSteps);
}

