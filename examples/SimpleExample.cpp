#include <iostream>
#include <string>
#include <cmath>

#include <Endas/Endas.hpp>
#include <Endas/Random/Random.hpp>
#include <Endas/DA/CovarianceOperator.hpp>
#include <Endas/DA/Algorithm/KalmanSmoother.hpp>


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
    endas::getRandomNumberGenerator().seed(1234);

    //-----------------------------------------------------------------------------------
    // Experiment setup
    //-----------------------------------------------------------------------------------

    // Number of data assimilation steps. The model has a period of 12 time steps.
    int nsteps = 12 * 10;

    // Observation error standard deviation
    double obsSigma = 0.4;

    // If greater than 1, observations will be assimilated on every n-th time step only
    int obsInterval = 3;


    // Model error covariance matrix
    DiagonalCovariance Q(makeArray({ 0.08, 0.01, 0.01 }).pow(2));

    // Initial (background) error covariance matrix. Purposefully quite large.
    DiagonalCovariance P0(makeArray({ 1.0, 1.0, 1.0 }).pow(2));

    // Observation error covariance matrix
    DiagonalCovariance R(makeArray({ obsSigma }).pow(2));

    // Observation operator. This is a 1x3 matrix representing a single observation of both
    // the first and second state variable. The third state variable is unobserved.
    MatrixObservationOperator H(makeMatrix(1, 3, { 1.0, 1.0, 0.0 }));

    // Initial system state
    Array x0 = makeArray({0, 0, 1});

    // Simple dynamic model for propagating the system state forward in time. The model is 
    // a rotation of the state around the first dimension by pi/6. The model is therefore 
    // periodic with 12 steps needed to complete one cycle.
    MatrixModel model(makeMatrix(3, 3, {
        1.0,  0.0,           0.0, 
        0.0,  cos(PI / 6.0), sin(PI / 6.0),
        0.0, -sin(PI / 6.0), cos(PI / 6.0)
    }));

    // Smoother lag. Zero disables smoothing (only filter solution is computed).
    int lag = endas::LAG_FIKS;

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
        tie(xtAll, zAll, obsTimes) = generateTestData(nsteps, x0, model, 1.0, H, Q, R, 
                                                      obsInterval);

        // Bad guess for the initial state
        Array x = makeArray({0.5, 0.5, 0.5});  
        
        // We will need P, Q and R as plain matrices for KalmanSmoother
        Matrix Pmat = P0.asMatrix();
        const Matrix& Qmat = Q.asMatrix();
        const Matrix& Rmat = R.asMatrix();
        const Matrix& Hmat = H.asMatrix();

        Array2d resultX(3, nsteps);
        Array2d resultSD(3, nsteps);
        resultX.col(0) = x;
        resultSD.col(0) = covError(Pmat);

        // This will be called on every KF/KS solution that becomes available
        kf.onResult([&](const Ref<const Array> x, const Ref<const Matrix> P, int k)
        {
            resultX.col(k) = x;
            resultSD.col(k) = covError(P);
        });
        
        //-----------------------------------------------------------------------------------
        // THE MAIN TIME-STEPPING LOOP 
        //-----------------------------------------------------------------------------------

        cout << "Running KS..." << endl;
        
        ENDAS_TIMER_BEGIN(ks);
        
        kf.beginSmoother(x, Pmat, 0);

        int obsIndex = 0;
        for (int k = 1; k != nsteps; k++)
        {
            // The Kalman Filter forecast step. The state `x` and error covariance `Pmat` are
            // updated in-place 
            kf.forecast(x, Pmat, Qmat, k-1, 1.0);

            // The update or analysis step
            kf.beginAnalysis(x, Pmat, k);

            // Do we have observations for this step? If yes, assimilate
            if (obsTimes[obsIndex] == k)
            {
                kf.assimilate(zAll.col(obsIndex++), Hmat, Rmat);
            }

            kf.endAnalysis();
        }

        kf.endSmoother();
        
        ENDAS_TIMER_END(ks);
        
        //-----------------------------------------------------------------------------------
        // Done, print statistics and generate plots 
        //-----------------------------------------------------------------------------------
        
        cout << "Kalman Filter/Smoother completed in " << ENDAS_TIMER_ELAPSED_MSEC(ks)
            << " milliseconds" << endl;

        Array resultRMSE = rmse(xtAll, resultX);
        cout << "Kalman Filter/Smoother mean RMSE: " << resultRMSE.mean() << endl;


#if ENDAS_PLOTTING_ENABLED

        auto xValues = rangeToVector(0, nsteps);

        map<string, string> lineStyle = { { "linewidth", "1"} };
        map<string, string> fillStyle = { { "alpha", "0.4"} };
        unordered_map<string, string> obsStyle  = { { "marker", "+"}, { "linewidth", "1"} };

        plt::figure_size(900, 900);
        
        for (int r = 0; r != 4; r++)
        {
            plt::subplot(4, 1, r+1);

            string ylabel;
            Array xt, xhat, xhat_sd;
            // First plot is X1+X2
            if (r == 0)
            {
                xt = xtAll.row(0) + xtAll.row(1);
                xhat = resultX.row(0) + resultX.row(1);
                xhat_sd = resultSD.row(0) + resultSD.row(1);
                ylabel = "X1+X2";
            }
            else
            {
                xt = xtAll.row(r-1);
                xhat = resultX.row(r-1);
                xhat_sd = resultSD.row(r-1);
                ylabel = "X" + to_string(r);
            }
            
            lineStyle["color"] = "black";
            plt::plot(xValues, toVector(xt), lineStyle);

            lineStyle["color"] = "green";
            fillStyle["color"] = "green";
            plt::plot(xValues, toVector(xhat), lineStyle);
            plt::fill_between(xValues, toVector(xhat - xhat_sd*1.96), 
                                       toVector(xhat + xhat_sd*1.96), 
                                       fillStyle);
            if (r == 0)
            {
                obsStyle["color"] = "purple";
                plt::scatter(obsTimes, toVector(zAll.row(0)), 13, obsStyle);
            }

            plt::xlabel("k");
            plt::ylabel(ylabel);
            plt::grid(true);
        }

        plt::show();

#endif

    }
    catch(const std::exception& e)
    {
        cerr << e.what() << endl;
    }


}