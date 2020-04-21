#include "Utils.hpp"

using namespace std;
using namespace endas;


tuple<Array2d, Array2d, vector<int>> 
endas::generateTestData(int nsteps, const Ref<const Array> x0, 
                        const GenericEvolutionModel& model, double dt,
                        const ObservationOperator& H, const CovarianceOperator& Q, 
                        const CovarianceOperator& R,
                        int obsInterval)
{
    int n = x0.size();
    int nobsSteps = (int)ceil((nsteps-1) / obsInterval);

    ENDAS_ASSERT(n > 0);
    ENDAS_ASSERT(nsteps > 0);

    // The 'true' system state evolution
    Array2d xtAll = Array2d::Zero(n, nsteps);
    xtAll.col(0) = x0;

    // Array of observations (we have nobs observations at each step 1..nsteps)
    Array2d zAll = Array2d::Zero(H.nobs(), nobsSteps);  

    // The observation times (steps)
    vector<int> obsSteps;

    Array x = x0;
    Array xnoise(n);
    Array z(H.nobs());
    Array znoise(H.nobs());

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
            H.apply(x, k, z);
            z += znoise;
            zAll.col(obsSteps.size()) = z;
            obsSteps.push_back(k);
        }
    }

    ENDAS_ASSERT(obsSteps.size() == zAll.cols());
    return make_tuple(xtAll, zAll, obsSteps);
}



std::tuple<Array2d, Array2d, double> 
endas::runKF(KalmanSmoother& kf, int nsteps, double dt, const Ref<const Array> x0, 
             const Ref<const Array2d> obs, const std::vector<int>& obsTimeSteps,
             const ObservationOperator& H, const CovarianceOperator& P0, 
             const CovarianceOperator& Q, const CovarianceOperator& R)
{
    int n = x0.size();
    Array x = x0;

    // We will need P, Q and R as plain matrices for KalmanSmoother
    Matrix Pmat, Qmat, Rmat, Hmat;
    P0.toMatrix(Pmat);
    Q.toMatrix(Qmat);
    R.toMatrix(Rmat);
    H.toMatrix(Hmat);

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
    ENDAS_TIMER_BEGIN(ks);

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
    
    ENDAS_TIMER_END(ks);
    double elapsed = ENDAS_TIMER_ELAPSED_SEC(ks);

    return make_tuple(move(resultX), move(resultSD), elapsed);
}



std::tuple<Array2d, Array2d, double> 
endas::runEnKF(EnsembleKalmanSmoother& kf, const GenericEvolutionModel& model,
               int nsteps, double dt, const Ref<const Array2d> E0, 
               const Ref<const Array2d> obs, const std::vector<int>& obsTimeSteps,
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
    ENDAS_TIMER_BEGIN(ks);

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
            kf.assimilate(obs.col(obsIndex++), H, R);
        }

        kf.endAnalysis();
    }

    kf.endSmoother();
    
    ENDAS_TIMER_END(ks);
    double elapsed = ENDAS_TIMER_ELAPSED_SEC(ks);

    return make_tuple(move(resultX), move(resultXSD), elapsed);
}

