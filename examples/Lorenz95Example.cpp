#include <iostream>
#include <string>
#include <cmath>
#include <chrono>

#include <Endas/Endas.hpp>
#include <Endas/Core/Profiling.hpp>
#include <Endas/Core/Ensemble.hpp>
#include <Endas/Random/Random.hpp>
#include <Endas/DA/CovarianceOperator.hpp>
#include <Endas/DA/Taper.hpp>
#include <Endas/DA/Algorithm/KalmanSmoother.hpp>
#include <Endas/DA/Algorithm/EnsembleKalmanSmoother.hpp>
#include <Endas/Domain/GenericStateSpace.hpp>

#include <EndasModels/Lorenz95.hpp>


#if ENDAS_PLOTTING_ENABLED
#   define WITHOUT_NUMPY
#   include "matplotlibcpp.hpp"
    namespace plt = matplotlibcpp;
#endif

#include "Utils.hpp"


using namespace std;
using namespace endas;


constexpr double PI = 3.141592653589793238462643383279502884;


int main(int argc, char *argv[])
{
    // Use pre-seeded RNG for deterministic output.
    // Please note that EnDAS should be configured with -DFORCE_DETERMINISTIC=ON to run 
    // deterministically
    endas::getRandomNumberGenerator().seed(1234);

    //-----------------------------------------------------------------------------------
    // Experiment setup
    //-----------------------------------------------------------------------------------

    // State size
    int n = 40;

    // Ensemble size
    int N = 20;

    // Forcing term of the Lorenz 96 model */
    int F = 8;

    // Number of data assimilation steps
    int nsteps = 1000;

    // Model integration time step.
    // Note: Too long integration step may cause difficulties with convergence, especially
    // for KalmanSmoother that relies on the linearization of the model. If this happens, 
    // decrease the integration step.
    double dt = 0.025 / 3.0;

    // Observation error standard deviation
    double obsSigma = 0.4;

    // If greater than 1, observations will be assimilated on every n-th time step only
    // Note: If this is set too high, filters will diverge. This is especially true of 
    // the KalmanSmoother. If this happens, lower this or decrease `dt`.
    int obsInterval = 6;


    double sigClim = 3.6414723;

    // Model error covariance matrix
    DiagonalCovariance Q(Array::Constant(n, 0.05 * sigClim).pow(2));

    // Initial (background) error covariance matrix
    DiagonalCovariance P0(Array::Constant(n, 0.5 * sigClim).pow(2));

    // Observation error covariance matrix
    DiagonalCovariance R(Array::Constant(n, 0.15 * sigClim).pow(2));

    // Observation operator. All state variables are observed
    int nobs = n;
    MatrixObservationOperator H(Matrix::Identity(n, n));

    // Observation 'coordinates' for use in localization. The index of the observed state 
    // variable will serve as the coordinate. We need 1 x nobs array
    Array2d obsCoords = makeSequence(0, nobs).transpose();


    // Initial system state with x21 perturbed a little
    Array x0 = Array::Constant(n, 8.0);
    x0[20] = 8.008;


    // Lorenz 96 nonlinear evolution model.
    Lorenz95Model model(n, F);

    // Smoother lag. Zero disables smoothing (only filter solution is computed)
    int lag = 10;

    // Kalman Filter we will be using
    endas::KalmanSmoother ks(model, lag);


    // Localization strategy for the ensemble smoothers:
    // Using generic one-dimensional state space where the index of each state variable is 
    // also its coordinate. Each state variable will be updated independently, i.e. in its 
    // own local domain.
    GenericStateSpace stateSpace(n);
    NoTaper taperFn(0);


    // Ensemble Kalman Filters we will be using
    endas::EnsembleKalmanSmoother enks(EnKS(), n, N, lag);
    enks.setCovInflationFactor(1.05);   // These are rough guesses, not tuned values.
    enks.localize(shared_ptr_wrap(stateSpace), shared_ptr_wrap(taperFn));

    // Ensemble Kalman Filters we will be using
    endas::EnsembleKalmanSmoother etks(ESTKS(), n, N, lag);
    etks.setCovInflationFactor(1.05);   // These are rough guesses, not tuned values.
    etks.localize(shared_ptr_wrap(stateSpace), shared_ptr_wrap(taperFn));

    //-----------------------------------------------------------------------------------
    // End of setup
    //-----------------------------------------------------------------------------------

    try
    {
        // Before we start, generate synthetic data that will serve as the "true" state and 
        // draw observations by applying the observation operator H and adding noise with
        // covariance R.
        cout << "Generating test data..." << endl;

        Array2d xtAll, obs;
        vector<int> obsTimes;

        tie(xtAll, obs, obsTimes) = generateTestData(nsteps, x0, model, dt, H, Q, R, 
                                                     obsInterval);

        // Use somewhat bad guess for x0 for data assimilation
        x0*= 1.5;


        // Generate initial ensemble from x0 and the background error covariance P0
        Array2d E0(n, N);
        generateEnsemble(x0, P0, E0);
        
        // Use the runXXX() utility functions to run the smoother time stepping loops

        cout << "Running KS..." << endl;
        Array2d ksX, ksErr;
        {
            ENDAS_PERF_SCOPE(KS);
            tie(ksX, ksErr) = runKF(ks, nsteps, dt, x0, obs, obsTimes, H, P0, Q, R);
        }

        cout << "Running EnKS..." << endl;
        Array2d enksX, enksErr;
        {
            ENDAS_PERF_SCOPE(ENKS);
            tie(enksX, enksErr) = runEnKF(enks, model, nsteps, dt, E0, obs, 
                                          obsCoords, obsTimes, H, Q, R);
        }

        cout << "Running ETKS..." << endl;
        Array2d etksX, etksErr;
        {
            ENDAS_PERF_SCOPE(ETKS);
            tie(etksX, etksErr) = runEnKF(etks, model, nsteps, dt, E0, obs, 
                                          obsCoords, obsTimes, H, Q, R);
        }


        // Done, print statistics and generate plots 

        Array ksRMSE = rmse(xtAll, ksX);
        Array enksRMSE = rmse(xtAll, enksX);
        Array etksRMSE = rmse(xtAll, etksX);

        cout << "KS mean RMSE  : " << ksRMSE.mean() << endl;
        cout << "EnKS mean RMSE: " << enksRMSE.mean() << endl;
        cout << "ETKS mean RMSE: " << etksRMSE.mean() << endl;

        cout << "Execution time summary:" << endl;
        profilingSummary(cout, 2);


#if ENDAS_PLOTTING_ENABLED

        // The plotted variable
        int X = 20;
        auto xValues = rangeToVector(0, nsteps);

        map<string, string> lineKeywords = { { "linewidth", "1"} };
        unordered_map<string, string> obsKeywords = { { "marker", "+"}, { "linewidth", "1"} };

        plt::figure_size(900, 600);

        plt::subplot(2, 1, 1);
        plt::title("Smoother performance on Lorenz96 problem");

        lineKeywords["color"] = "black";
        lineKeywords["label"] = "Truth";
        plt::plot(xValues, toVector(xtAll.row(20)), lineKeywords);

        obsKeywords["color"] = "purple";
        obsKeywords["label"] = "Observation";
        plt::scatter(obsTimes, toVector(obs.row(X)), 13, obsKeywords);

        lineKeywords["color"] = "green";
        lineKeywords["label"] = "KS";
        plt::plot(xValues, toVector(ksX.row(X)), lineKeywords);
        
        lineKeywords["color"] = "red";
        lineKeywords["label"] = "EnKS";
        plt::plot(xValues, toVector(enksX.row(X)), lineKeywords);
        
        lineKeywords["color"] = "blue";
        lineKeywords["label"] = "ESTKS";
        plt::plot(xValues, toVector(etksX.row(X)), lineKeywords);
        
        plt::xlabel("k");
        plt::ylabel("X" + to_string(X));
        plt::grid(true);
        plt::legend();
        

        plt::subplot(2, 1, 2);
        lineKeywords["color"] = "green";
        lineKeywords["label"] = "KS";
        plt::plot(xValues, toVector(ksRMSE), lineKeywords);

        lineKeywords["color"] = "red";
        lineKeywords["label"] = "EnKS";
        plt::plot(xValues, toVector(enksRMSE), lineKeywords);

        lineKeywords["color"] = "blue";
        lineKeywords["label"] = "ESTKS";
        plt::plot(xValues, toVector(etksRMSE), lineKeywords);

        plt::xlabel("k");
        plt::ylabel("RMSE");
        plt::grid(true);
        plt::legend();

        plt::show();

#endif

    }
    catch(const std::exception& e)
    {
        cerr << e.what() << endl;
    }


}