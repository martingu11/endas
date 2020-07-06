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
#include <Endas/DA/GridDomain.hpp>
#include <Endas/DA/GridDomainPartitioning.hpp>
#include <Endas/DA/SimpleObservationManager.hpp>

#include <Endas/IO/ArrayIO.hpp>
#include <Endas/Utils/EnsembleSampler.hpp>
#include <Endas/Spatial/Variogram.hpp>
#include <Endas/Random/GaussianRandomField.hpp>
#include <Endas/Random/MultivariateRandomNormal.hpp>

#include <EndasModels/QG.hpp>


#if ENDAS_PLOTTING_ENABLED
#   include "matplotlibcpp.hpp"
    namespace plt = matplotlibcpp;
#endif

#include "Utils.hpp"

using namespace std;
using namespace endas;


/*
 * Custom observation operator used in the QG example.
 * The operator systematically selects observations such that state variables x = o + m*i are 
 * directly observed, where o is an offset, m is the "spacing" and i = 0,1,2,... The operator 
 * implements subset() so that it can be used with localized algorithms. We do not need to 
 * implement toDenseMatrix().
 */
class QGObservationOperator : public ObservationOperator
{
public:

    QGObservationOperator(const GriddedDomain& domain, index_t m, index_t offset)
    : mDomain(domain),
      mStateSize(domain.size()), mM(m),
      mOffset(offset % m)
    { }

    virtual index_t nobs() const 
    { 
        return (index_t)ceil((mStateSize - mOffset) / (double)mM); 
    }

    virtual index_t nstate() const { return mStateSize; }
    virtual bool isLinear() const { return true; }

    virtual void apply(const Ref<const Array2d> x, Ref<Array2d> out) const
    {
        ENDAS_ASSERT(x.rows() == mStateSize);
        ENDAS_ASSERT(out.rows() == nobs());

        // Select every m-th variable, starting with an offset
        index_t i = 0;
        for (index_t j = mOffset; j < mStateSize; j+= mM) 
            out.row(i++) = x.row(j);
    }

    virtual shared_ptr<const ObservationOperator> subset(const IndexArray& indices) const
    {
        // We need to construct an observation operator that will represent a subset of the 
        // observation space. We can use CustomObservationOperator for this. Please note that
        // the lambda capture must be by value to avoid segmentation faults!

        index_t nobs = indices.size();

        return make_shared<CustomObservationOperator>(indices.size(), mStateSize, true, 
        [=](const Ref<const Array2d> x, Ref<Array2d> out)
        {
            ENDAS_ASSERT(out.rows() == nobs);
            for (index_t i = 0; i != nobs; i++) 
                out.row(i) = x.row(indices[i]*mM + mOffset);
        });
    }

    
    // Returns coordinates of the observations. These will be needed for localization to work
    // via the SimpleObservationManager. 
    Array2d obsCoords() const
    {
        IndexArray obsIndices;
        for (index_t i = mOffset; i < mStateSize; i+= mM) obsIndices.push_back(i);

        Array2d out(mDomain.coordDim(), obsIndices.size());
        mDomain.getCoords(obsIndices, out);
        return move(out);
    }

private:
    index_t mOffset, mM, mStateSize;
    const GriddedDomain& mDomain;
};


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
    int numSpinupSteps = 10000;

    // The number of model steps used for generation of the initial ensemble
    int numEnsembleInitSteps = 1000;

    // Number of data assimilation steps
    int numAssimSteps = 20;

    // 1.5 -layer quasi-geostrophic model with internal time-stepping `modeldt`.
    QGModel model(N, modeldt);
    int nx = model.sizex();
    int ny = model.sizey();
    int n = ny*nx;    

    // Observation error variance
    double obsVariance = 4.0;

    // Initial mean for the data assimilation run
    double initMeanDA = 0;

    // We aim at having roughly 300 observations of the state space at every assimilation step
    int numObs = 300;

    // If greater than 1, observations will be assimilated on every n-th time step only
    // Note: If this is set too high, filters may diverge. 
    int obsInterval = 4;


    // Zero model error covariance matrix (i.e. perfect model scenario)
    ZeroCovariance Q(n);

   
    // State space is a two-dimensional grid in Euclidean space. The size of the grid is
    // ny*nx and it occupies space (0, 0, nx, ny). Therefore, each cell is 1x1 "spatial unit"
    // large.
    GridDomain stateSpace(makeShape(ny, nx), make_shared<EuclideanCS>(2), 
                          makeBox<AABox>(0, 0, ny, nx));

    // We will partition the state space grid into local blocks for localized analysis.
    // This is necessary, without localization the filter will diverge after few dozens steps.
    // We use blocks of 3x3 cells as local analysis domains. We will use Gaspari-Cohn tapering
    // function to limit the effect of observations as distance increases. The de-corellation 
    // distance will be set to 15 (cells).

    GridDomainPartitioning statePartitioner(shared_ptr_wrap(stateSpace), 3);
    SphericalTaper taperFn(15);

    // Smoother lag. Zero disables smoothing (only filter solution is computed).
    int lag = 10;

    // Ensemble Kalman Filters we will be using
    endas::EnsembleKalmanSmoother enks(ESTKS(), n, N, lag);
    enks.setCovInflationFactor(1.05);   // These are rough guesses, not tuned values.
    enks.localize(shared_ptr_wrap(statePartitioner), shared_ptr_wrap(taperFn));

    // File storing the initial ensemble in NumPy array format
    string xt0FilePath = "./qgexample_xt0.npy";
    string SFilePath = "./qgexample_S.npy";
    string UFilePath = "./qgexample_U.npy";


    //-----------------------------------------------------------------------------------
    // End of setup
    //-----------------------------------------------------------------------------------

    try
    {

        // Generate the initial state of the system using a long model run. Because this 
        // does take time, we will store the state to disk for reuse
        SecondOrderExactSample sampler;
        Array2d xt0 = Array::Zero(n);

        if (!exists(xt0FilePath) || !exists(SFilePath) || !exists(UFilePath))
        {
            cout << "Performing long model run to generate initial state and background ";
            cout << "covariance data. This may take a while..." << endl;
            
            Array2d states = Array2d::Zero(n, numEnsembleInitSteps);

            QGModel model(1, modeldt);
            for (int i = 0; i != numSpinupSteps; i++)
            {
                model(xt0, i, modeldt);
            }

            for (int i = 0; i != numEnsembleInitSteps; i++)
            {
                model(xt0, numSpinupSteps+i, modeldt);
                states.col(i) = xt0;
            }

            cout << "Model run completed, computing EOFs..." << endl;
            sampler.initFromStates(states, true, 0.001);
            cout << "Have " << sampler.numEOFs() << " EOFs..." << endl;

            saveArrayAsNpy(xt0, xt0FilePath);     
            saveArrayAsNpy(sampler.getS(), SFilePath);       
            saveArrayAsNpy(sampler.getU(), UFilePath);       
            cout << "Initial state and EOFs saved" << endl;
        }
        else 
        {
            cout << "Using initial state and covariance data from files" << endl;
            xt0 = loadArrayFromNpy(xt0FilePath);
            auto S = loadArrayFromNpy(SFilePath);
            auto U = loadArrayFromNpy(UFilePath);
            sampler.initFromEOFs(move(S), move(U));
        }

        // Generate the true state trajectory and observations first. To do that,
        // we will be using different observation and covariance operator for every time
        // step with a random offset to the observed states. Because we will need the 
        // operators also later in the assimilation round, we store the offsets.
        
        cout << "Generating true state trajectory and observations..." << endl;

        int HOperatorSpacing = (int)(n / (double)numObs);
        unordered_map<int, int> HOperatorOffsets;

        Array2d xtAll;
        vector<Array> obsAll;
        tie(xtAll, obsAll) = generateExampleData(
            numAssimSteps, xt0, 
            QGModel(1, modeldt), modeldt, 
            Q,
            [&](int k)
            {
                int offset = rand() % HOperatorSpacing;
                HOperatorOffsets[k] = offset;

                auto H = make_shared<QGObservationOperator>(
                    stateSpace, HOperatorSpacing, offset);

                auto R = make_shared<DiagonalCovariance>(H->nobs(), obsVariance);
                return ObservationOpAndCov { H, R};
            },
            0, obsInterval);
        

        cout << "Generating initial ensemble..." << endl;

        Array2d E0(n, N);
        sampler.samplePerturbations(E0);
        saveArrayAsNpy(E0, "./qgexample_E0.npy");


        //-----------------------------------------------------------------------------------
        // Ensemble Kalman Smoother time-stepping loop
        //-----------------------------------------------------------------------------------

        cout << "Running EnKS..." << endl;

        Array2d valuesEnKS(n, numAssimSteps);
        Array2d errorsEnKS(n, numAssimSteps);
        valuesEnKS.col(0) = ensembleMean(E0);
        errorsEnKS.col(0) = ensembleError(E0);

        // This will be called on every KF/KS solution that becomes available
        enks.onResult([&](const Ref<const Array2d> E, int k)
        {
            valuesEnKS.col(k) = ensembleMean(E);
            errorsEnKS.col(k) = ensembleError(E);
        });
        
        // The main time-stepping loop starts here
        Array2d E = E0;
        enks.beginSmoother(E, 0);

        for (int k = 1; k != numAssimSteps; k++)
        {
            // Propagate the ensemble from time step k-1 to k
            model(E, k, modeldt);

            // The update or analysis step
            enks.beginAnalysis(E, k);

            // Do we have observations for this step? If yes, assimilate
            if (obsAll[k].size() > 0)
            {
                cout << "Assimilating " << obsAll[k].size() << " observations at time step " << k << "..." << endl;

                // We need to construct the observation operator and error covariance 
                // identical to those used in generateExampleData()
                auto H = make_shared<QGObservationOperator>(
                    stateSpace, HOperatorSpacing, HOperatorOffsets[k]);

                auto R = make_shared<DiagonalCovariance>(H->nobs(), obsVariance);
                ENDAS_ASSERT(H->nobs() == obsAll[k].size());

                SimpleObservationManager omgr(obsAll[k], H->obsCoords(), H, R);
                enks.assimilate(omgr);
            }

            enks.endAnalysis();
        }

        enks.endSmoother();

        // Done with smoothing, plot some nice figures

#if ENDAS_PLOTTING_ENABLED

        auto xValues = rangeToVector(0, numAssimSteps);

        // Plotting styles
        unordered_map<string, string> obsstyle = { {"marker","+"}, {"linewidth","1"} };
        map<string, string> xstyle = { {"vmin","-40"}, {"vmax","40"}, {"cmap","Spectral"} };
        map<string, string> errorstyle = { {"vmin","0"}, {"vmax","4"}, {"cmap","YlGnBu"} };

        // Which tim steps we will plot? We want few steps over the assimilation time window, 
        // those for which we had observations and the step k=0.
        vector<int> Ks;        
        int numPlots = 4; // Print plots for 4 time steps

        double K = numAssimSteps / (double)numPlots;  
        for (int i = 0; i != numPlots; i++)
        {
            int k = (int)((K * (double)i) / obsInterval) * obsInterval;
            Ks.push_back(k);
        }

        int nrows = 3;
        int ncols = Ks.size();
        int col = 1;

        plt::figure_size(380*ncols, 300*nrows);
        PyObject* fig;

        for (int k : Ks)
        {
            stringstream title;
            title << "k = " << k;

            // Plot the true field and locations of the generated observations
            plt::subplot(nrows, ncols, col);
            plt::title(title.str());
            plt::imshow(xtAll.col(k).data(), ny, nx, 1, xstyle, &fig);
            if (k == 0) plt::ylabel("True state and observations");
            if (k == Ks.back()) plt::colorbar(fig);

            if (obsAll[k].size() > 0)
            {
                QGObservationOperator H(stateSpace, HOperatorSpacing, HOperatorOffsets[k]);
                auto coords = H.obsCoords();
                plt::scatter(toVector(coords.row(1)), toVector(coords.row(0)), 1, obsstyle);
            }

            // Plot the EnkF/EnkS fields
            plt::subplot(nrows, ncols, col + ncols*1);
            plt::imshow(valuesEnKS.col(k).data(), ny, nx, 1, xstyle, &fig);
            if (k == 0) plt::ylabel("Estimate");
            if (k == Ks.back()) plt::colorbar(fig);

            // Plot the ensemble spread
            plt::subplot(nrows, ncols, col + ncols*2);
            plt::imshow(errorsEnKS.col(k).data(), ny, nx, 1, errorstyle, &fig);
            if (k == 0) plt::ylabel("Ensemble spread (1 STD)");
            if (k == Ks.back()) plt::colorbar(fig);

            ++col;
        }

        plt::tight_layout();
        plt::show();

#endif

    }
    catch(const std::exception& e)
    {
        cerr << e.what() << endl;
    }


}