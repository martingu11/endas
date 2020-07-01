#include "Utils.hpp"
#include <Endas/Core/Ensemble.hpp>
#include <Endas/DA/ObservationOperator.hpp>
#include <Endas/DA/CovarianceOperator.hpp>
#include <Endas/Domain/SimpleObservationManager.hpp>

using namespace std;
using namespace endas;


tuple<Array2d, vector<Array>> 
endas::generateExampleData(int numSteps, const Ref<const Array> x0, 
                           const EvolutionModel& model, double dt,
                           const ObservationOperator& H, const CovarianceOperator& Q, 
                           const CovarianceOperator& R,
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
        Array znoise(H.nobs());
        R.randomMultivariateNormal(znoise);

        if (k % obsInterval == 0)
        {
            Array z(H.nobs());
            H.apply(x, z);
            z += znoise;
            zAll.push_back(move(z));
        }
        else 
            zAll.push_back(Array());
    }

    //ENDAS_ASSERT(obsSteps.size() == zAll.cols());
    return make_tuple(xtAll, zAll);
}



/*std::tuple<Array2d, Array2d> 
endas::runKF(KalmanSmoother& kf, int nsteps, double dt, const Ref<const Array> x0, 
             const Ref<const Array2d> obs, const std::vector<int>& obsTimeSteps,
             const ObservationOperator& H, const CovarianceOperator& P0, 
             const CovarianceOperator& Q, const CovarianceOperator& R)
{
    int n = x0.size();
    Array x = x0;

    // We will need P, Q and R as plain matrices for KalmanSmoother
    Matrix Pmat = P0.toDenseMatrix();
    const Matrix& Qmat = Q.toDenseMatrix();
    const Matrix& Rmat = R.toDenseMatrix();
    const Matrix& Hmat = H.toDenseMatrix();

    Array2d resultX(n, nsteps);
    Array2d resultSD(n, nsteps);
    resultX.col(0) = x;
    resultSD.col(0) = covError(Pmat);

    // This will be called on every KF/KS solution that becomes available
    kf.onResult([&](const Ref<const Array> x, const Ref<const Matrix> P, int k)
    {
        resultX.col(k) = x;
        resultSD.col(k) = covError(P);
    });
    
    // The main time stepping loop
    kf.beginSmoother(x, Pmat, 0);

    int obsIndex = 0;
    for (int k = 1; k != nsteps; k++)
    {
        // The Kalman Filter forecast step. The state `x` and error covariance `Pmat` are
        // updated in-place 
        kf.forecast(x, Pmat, Qmat, k-1, dt);

        // The update or analysis step
        kf.beginAnalysis(x, Pmat, k);

        // Do we have observations for this step? If yes, assimilate
        if (obsTimeSteps[obsIndex] == k)
        {
            kf.assimilate(obs.col(obsIndex++), Hmat, Rmat);
        }

        kf.endAnalysis();
    }

    kf.endSmoother();
    
    return make_tuple(move(resultX), move(resultSD));
}



std::tuple<Array2d, Array2d> 
endas::runEnKF(EnsembleKalmanSmoother& kf, const EvolutionModel& model,
               int nsteps, double dt, const Ref<const Array2d> E0, 
               const Ref<const Array2d> obs, const Ref<const Array2d> obsCoords, 
               const std::vector<int>& obsTimeSteps,
               const ObservationOperator& H, const CovarianceOperator& Q, 
               const CovarianceOperator& R)
{
    int n = E0.rows();
    int N = E0.cols();

    Array2d E = E0;

    // We will only store ensemble mean and error(SD)
    Array2d resultX(n, nsteps);
    Array2d resultXSD(n, nsteps);
    
    resultX.col(0) = ensembleMean(E0);
    resultXSD.col(0) = ensembleError(E0);

    // This will be called on every KF/KS solution that becomes available
    kf.onResult([&](const Ref<const Array2d> E, int k)
    {
        resultX.col(k) = ensembleMean(E);
        resultXSD.col(k) = ensembleError(E);
    });
    
    // The main time stepping loop
    kf.beginSmoother(E, 0);

    int obsIndex = 0;
    for (int k = 1; k != nsteps; k++)
    {
        // The Ensemble Kalman Filter forecast step. 
        ensembleForecast(E, model, Q, k, dt);

        // The update or analysis step
        kf.beginAnalysis(E, k);

        // Do we have observations for this step? If yes, assimilate
        if (obsTimeSteps[obsIndex] == k)
        {
            Ref<const Array> z = obs.col(obsIndex++);
            SimpleObservationManager omgr(z, obsCoords, shared_ptr_wrap(H), shared_ptr_wrap(R));
            kf.assimilate(omgr);
        }

        kf.endAnalysis();
    }

    kf.endSmoother();
    
    return make_tuple(move(resultX), move(resultXSD));
}*/

