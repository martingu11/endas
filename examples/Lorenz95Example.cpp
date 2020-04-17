#include <iostream>
#include <string>
#include <cmath>
#include <chrono>

#include <Endas/Endas.hpp>
#include <Endas/Random/Random.hpp>
#include <Endas/Algorithm/KalmanSmoother.hpp>
#include <Endas/Error/DiagonalCovariance.hpp>
#include <Endas/Model/Lorenz95.hpp>


#if ENDAS_PLOTTING_ENABLED
#   define WITHOUT_NUMPY
#   include "matplotlibcpp.hpp"
    namespace plt = matplotlibcpp;
#endif

#include "Utils.hpp"


using namespace std;
using namespace endas;

typedef chrono::steady_clock perfclock_t;

constexpr double PI = 3.141592653589793238462643383279502884;


int main(int argc, char *argv[])
{
    endas::getRngEngine().seed(1234);

    //-----------------------------------------------------------------------------------
    // Experiment setup
    //-----------------------------------------------------------------------------------

    int n = 40;

    int F = 8;

    // Number of data assimilation steps
    int nsteps = 1000;

    // Model integration time step (equivalent to 30min)
    double dt = 0.025 / 6.0;

    // Observation error standard deviation
    double obsSigma = 0.4;

    // If greater than 1, observations will be assimilated on every n-th time step only
    int obsInterval = 6;


    double sigClim = 3.6414723;

    // Model error covariance matrix
    DiagonalCovariance Q(Array::Constant(n, 0.05 * sigClim).pow(2));

    // Initial (background) error covariance matrix
    DiagonalCovariance P0(Array::Constant(n, 0.5 * sigClim).pow(2));

    // Observation error covariance matrix
    DiagonalCovariance R(Array::Constant(n, 0.15 * sigClim).pow(2));

    // Observation operator. We will observe the last 3 variables in every 5 state vector 
    // variables, i.e. 24 out of 40 in this case
    int nobs = 40;
    Matrix H = Matrix::Identity(n, n);

    // Initial system state with x21 perturbed a little
    Array x0 = Array::Constant(n, 8.0);
    x0[20] = 8.008;


    // Simple dynamic model for propagating the system state forward in time. The model is 
    // a rotation of the state around the first dimension by pi/6. The model is therefore 
    // periodic with 12 steps needed to complete one cycle.
    Lorenz95Model model(n, F);

    // Smoother lag. Zero disables smoothing (only filter solution is computed)
    int lag = 1;

    // Kalman Filter we will be using
    endas::KalmanSmoother kf(model, lag);

    //-----------------------------------------------------------------------------------
    // End of setup
    //-----------------------------------------------------------------------------------

    try
    {
        // Before we start, generate synthetic data that will serve as the "true" state and 
        // draw observations by applying the observation operator H and adding noise with
        // covariance R.
        cout << "Generating test data..." << endl;

        Array2d xtAll, zAll;
        vector<int> obsTimes;

        tie(xtAll, zAll, obsTimes) = generateTestData(nsteps, x0, model, dt, H, Q, R, obsInterval);

        // Bad guess for the initial state
        Array x = x0;
        
        // We will need P, Q and R as plain matrices for KalmanSmoother
        Matrix Pmat, Qmat, Rmat;
        P0.toMatrix(Pmat);
        Q.toMatrix(Qmat);
        R.toMatrix(Rmat);

        Array2d resultX(n, nsteps);
        Array2d resultSD(n, nsteps);
        resultX.col(0) = x;
        resultSD.col(0) = cov2error(Pmat);

        // This will be called on every KF/KS solution that becomes available
        kf.onResult([&](const Ref<const Array> x, const Ref<const Matrix> P, int k)
        {
            resultX.col(k) = x;
            resultSD.col(k) = cov2error(P);

            //cout << "onresult " << k << " (" << x.transpose() << ")" << endl;
        });
        
        //-----------------------------------------------------------------------------------
        // THE MAIN TIME-STEPPING LOOP 
        //-----------------------------------------------------------------------------------

        cout << "Running KF..." << endl;

        endas::getRngEngine().seed(1234);
        auto timerBegin = perfclock_t::now();

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
            if (obsTimes[obsIndex] == k)
            {
                kf.assimilate(zAll.col(obsIndex++), H, Rmat);
            }

            kf.endAnalysis();
        }

        kf.endSmoother();
        auto timerEnd = perfclock_t::now();

        //-----------------------------------------------------------------------------------
        // Done, print statistics and generate plots 
        //-----------------------------------------------------------------------------------

        cout << "Kalman Filter/Smoother completed in " 
            << chrono::duration_cast<chrono::milliseconds>(timerEnd - timerBegin).count() / 1000.0 
            << " seconds" << endl;


        Array resultRMSE = rmse(xtAll, resultX);
        cout << "Kalman Filter/Smoother mean RMSE: " << resultRMSE.mean() << endl;


#if ENDAS_PLOTTING_ENABLED

        // The plotted variable
        int X = 20;

        auto xValues = rangeToVector(0, nsteps);

        map<string, string> lineStyle = { { "linewidth", "1"} };
        map<string, string> fillStyle = { { "alpha", "0.4"} };
        unordered_map<string, string> obsStyle  = { { "marker", "+"}, { "linewidth", "1"} };

        plt::figure_size(900, 600);
        
        lineStyle["color"] = "black";
        plt::plot(xValues, toVector(xtAll.row(20)), lineStyle);

        lineStyle["color"] = "green";
        fillStyle["color"] = "green";
        plt::plot(xValues, toVector(resultX.row(X)), lineStyle);
        plt::fill_between(xValues, toVector(resultX.row(X) - resultSD.row(X)*1.96), 
                                    toVector(resultX.row(X) + resultSD.row(X)*1.96), 
                                    fillStyle);
        
        obsStyle["color"] = "purple";
        plt::scatter(obsTimes, toVector(zAll.row(X)), 13, obsStyle);

        /*plt::xlabel("k");
        plt::ylabel("X" + to_string(r+1));
        plt::grid(true);
        

        plt::subplot(4, 1, 4);
        lineStyle["color"] = "black";
        plt::plot(xValues, toVector(resultRMSE), lineStyle);
        plt::xlabel("k");
        plt::ylabel("RMSE");
        plt::ylim(0.0, obsSigma);
        plt::grid(true);*/

        plt::show();

#endif

    }
    catch(const std::exception& e)
    {
        cerr << e.what() << endl;
    }


}