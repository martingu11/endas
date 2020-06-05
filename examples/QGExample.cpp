/**
 * @file QGExample.cpp
 * 
 * Data assimilation example using the 1.5 -layer Quasi-Geostrophic model.
 * For more information about the model and the data assimilation set up, see 
 * 
 * SAKOV, P. and OKE, P.R. (2008), A deterministic formulation of the ensemble Kalman filter: 
 * an alternative to ensemble square root filters. Tellus A, 60: 361-371. 
 * doi:10.1111/j.1600-0870.2007.00299.x
 */

#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <sys/stat.h>

#include <Endas/Endas.hpp>
#include <Endas/Core/Profiling.hpp>
#include <Endas/Core/Ensemble.hpp>
#include <Endas/Random/Random.hpp>
#include <Endas/DA/CovarianceOperator.hpp>
#include <Endas/DA/Taper.hpp>
#include <Endas/DA/Algorithm/KalmanSmoother.hpp>
#include <Endas/DA/Algorithm/EnsembleKalmanSmoother.hpp>
#include <Endas/Domain/GridStateSpace.hpp>

#include <EndasModels/QG.hpp>


#if ENDAS_PLOTTING_ENABLED
#   include "matplotlibcpp.hpp"
    namespace plt = matplotlibcpp;
#endif

#include "Utils.hpp"

using namespace std;
using namespace endas;


// Checks whether file exists
inline bool exists(const std::string& name) 
{
    struct stat buffer;   
    return (stat (name.c_str(), &buffer) == 0); 
}


int main(int argc, char *argv[])
{
    // Use pre-seeded RNG for deterministic output.
    // Please note that EnDAS should be configured with -DFORCE_DETERMINISTIC=ON to run 
    // deterministically
    endas::getRandomNumberGenerator().seed(1234);

    //-----------------------------------------------------------------------------------
    // Experiment setup
    //-----------------------------------------------------------------------------------

    // Ensemble size
    int N = 25;

    // Internal QG model integration time step
    double modeldt = 1.25;

    // Data assimilation is performed every 4th time step
    double assimInterval = 4;

    // Number of model integration steps for generating the initial system state
    int numInitSteps = 10000;

    // Number of data assimilation steps
    int numAssimSteps = 300;

    // Observation error variance
    double obsVariance = 4.0;

    // If greater than 1, observations will be assimilated on every n-th time step only
    // Note: If this is set too high, filters will diverge. This is especially true of 
    // the KalmanSmoother. If this happens, lower this or decrease `dt`.
    //int obsInterval = 6;

    // 1.5 -layer quasi-geostrophic model with internal time-stepping `dt`.
    QGModel model(N, modeldt);

    // Smoother lag. Zero disables smoothing (only filter solution is computed)
    int lag = 10;

    // Localization strategy for the ensemble smoothers:
    // Using generic one-dimensional state space where the index of each state variable is also its 
    // coordinate. Each state variable will be updated independently, i.e. in its own local domain.
    GridStateSpace stateSpace();
    NoTaper taperFn(0);

    int nx = model.sizex();
    int ny = model.sizey();
    int n = ny*nx;

    // Ensemble Kalman Filters we will be using
    endas::EnsembleKalmanSmoother enks(EnKS(), n, N, lag);
    enks.setCovInflationFactor(1.05);   // These are rough guesses, not tuned values.
    //enks.localize(shared_ptr_wrap(stateSpace), shared_ptr_wrap(taperFn));

    // Ensemble Kalman Filters we will be using
    endas::EnsembleKalmanSmoother etks(ESTKS(), n, N, lag);
    etks.setCovInflationFactor(1.05);   // These are rough guesses, not tuned values.
    //etks.localize(shared_ptr_wrap(stateSpace), shared_ptr_wrap(taperFn));

    // File storing the initial system state in NumPy array format
    string initialStateFile = "./qgexample_x0.npy";


    //-----------------------------------------------------------------------------------
    // End of setup
    //-----------------------------------------------------------------------------------

    try
    {
        // Generate the initial state of the system using a long model run. Because this does
        // take time, we will store the state to disk for reuse
        Array x0;
        if (!exists(initialStateFile))
        {
            cout << "Generating initial system state. This may take a while..." << endl;
            x0 = Array::Zero(n);

            QGModel model(N, modeldt);
            model(x0, 0, modeldt * numInitSteps);
            saveAsNpy(x0, initialStateFile);       
            cout << "Initial system state saved to " << initialStateFile << endl;
        }
        else 
        {
            cout << "Using initial state data from " << initialStateFile << endl;
            x0 = loadFromNpy(initialStateFile);
            ENDAS_ASSERT(x0.size() == n);
        }


        // Model error covariance matrix of a perfect model
        auto Q = ZeroCovariance(n);



        // Before we start, generate synthetic data that will serve as the "true" state and 
        // draw observations by applying the observation operator H and adding noise with
        // covariance R.
        //cout << "Generating test data for " << numDAsteps << " steps..." << endl;

        //Array2d xtAll, obs;
        //vector<int> obsTimes;


        // Observation error covariance matrix
        //DiagonalCovariance R(Array::Constant(n, 4.0));


        

        // Observation operator. All state variables are observed
        //int nobs = n;
        //CustomObservationOperator H(n, n, true, [](const Ref<const Array2d> x, int k, Ref<Array2d> out)
        //{
        //    out = x; 
        //});

        // Observation 'coordinates' for use in localization. The index of the observed state variable
        // will serve as the coordinate. We need 1 x nobs array
        //Array2d obsCoords = makeSequence(0, nobs).transpose();

        

#if 0
//#if ENDAS_PLOTTING_ENABLED

        // The plotted variable
        int X = 20;
        auto xValues = rangeToVector(0, nsteps);

        map<string, string> lineKeywords = { { "linewidth", "1"} };
        unordered_map<string, string> obsKeywords = { { "marker", "+"}, { "linewidth", "1"} };

        plt::figure_size(900, 600);

        plt::subplot(2, 1, 1);
        plt::title("Smoother performance on Lorenz96 problem");

        Array psi(n);
        //model.calc_psi(xtAll.col(nsteps-1), psi); 
        //plt::imshow(psi.data(), ny, nx, 1);
        cout << xtAll.col(nsteps-1).mean() << endl;

        plt::imshow(xtAll.col(nsteps-1).data(), ny, nx, 1);


        /*plt::subplot(2, 1, 2);
        model.calc_psi(xtAll.col(nsteps-2), psi); 
        plt::imshow(psi.data(), ny, nx, 1);*/


        /*lineKeywords["color"] = "black";
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
        plt::legend();*/

        plt::show();

#endif

    }
    catch(const std::exception& e)
    {
        cerr << e.what() << endl;
    }


}