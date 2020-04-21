#include <Endas/Algorithm/EnsembleKalmanSmoother.hpp>
#include <Endas/Caching/MemoryArrayCache.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"

#include <vector>
#include <deque>

using namespace std;
using namespace endas;

using handle_t = ArrayCache::handle_t;

typedef SequentialEnsembleSmoother::OnResultFn OnResultFn;


EnKSVariant::~EnKSVariant() { }

void EnKSVariant::init(int n, int N)
{ }

void EnKSVariant::applyCovInflation(Ref<Array2d> E, double factor, int k) const
{
    inflateInPlace(E, factor);
}

void EnKSVariant::processGlobalEnsemble(const Ref<const Array2d> Ag, const ObservationOperator& H,
                                        int k, vector<Array2d>& dataOut) const
{ }

/* Smoother data for one time step. */
struct ENKSData
{
    int k;
    handle_t E;
};

struct EnsembleKalmanSmoother::Data
{
    unique_ptr<EnKSVariant> variant;

    int n, N;
    int lag;
    double covInflation;

    // Localiaztion
    int numDomains;

    // Data for the current update
    bool updateActive;
    int upK;
    unique_ptr<Ref<Array2d>> upE; // Have to use pointer to Ref<> until we have Ptr<>
    Matrix upX; // Ensemble transform matrix
    bool upHaveX; // Is upX already holding a transform?

    
    // Smoother data
    double forgettingFactor;
    shared_ptr<ArrayCache> cache;
    deque<ENKSData> ksSteps; // Data for all smoother steps

    Data() 
    : lag(0), covInflation(1.0), forgettingFactor(1.0), 
      updateActive(false), numDomains(0)
    { }

    Data(const Data&) = delete;

    void pushSmootherStep(int k, handle_t E)
    {
        ksSteps.push_back({ k, E });
    }

    //void assimilate(Ref<Array2d> E, void* Egdata, const Ref<const Array> z, 
    //                const ObservationOperator& H, const CovarianceOperator& R);

    void laggedSmoother(OnResultFn& onresult, bool finishing);

};


EnsembleKalmanSmoother::EnsembleKalmanSmoother(const EnKSVariant& variant,  int n, int N,
                                               int lag, shared_ptr<ArrayCache> cache)
: mData(make_unique<Data>())
{
    ENDAS_ASSERT(lag >= 0);

    mData->n = n;
    mData->N = N;
    mData->variant = variant.clone();
    mData->lag = lag;
    mData->cache = (cache)? cache : make_shared<MemoryArrayCache>();

    mData->variant->init(n, N);
}


EnsembleKalmanSmoother::~EnsembleKalmanSmoother()
{ }


void EnsembleKalmanSmoother::setCovInflationFactor(double factor)
{
    ENDAS_ASSERT(factor >= 1.0);
    mData->covInflation = factor;
}

void setSmootherForgettingFactor(double factor);

void EnsembleKalmanSmoother::setSmootherForgettingFactor(double factor)
{
    ENDAS_ASSERT(factor <= 1.0);
    mData->forgettingFactor = factor;
}


void EnsembleKalmanSmoother::beginSmoother(const Ref<const Array2d> E0, int k0)
{
    ENDAS_ASSERT(!mData->updateActive);
    ENDAS_ASSERT(E0.rows() == mData->n);
    ENDAS_ASSERT(E0.cols() == mData->N);

    if (mData->lag == 0) return;

    mData->ksSteps.clear();
    mData->upE.reset();
    
    auto E0handle = mData->cache->put(E0);
    mData->pushSmootherStep(k0, E0handle);
}



void EnsembleKalmanSmoother::beginAnalysis(Ref<Array2d> E, int k)
{
    ENDAS_ASSERT(!mData->updateActive);
    ENDAS_ASSERT(E.rows() == mData->n);
    ENDAS_ASSERT(E.cols() == mData->N);

    int N = E.cols();

    mData->upK = k;
    mData->upE = make_unique<Ref<Array2d>>(E);
    mData->upX.resize(N, N);
    mData->upHaveX = false;

    // Multiplicative covariance inflation 
    if (mData->covInflation != 1.0)
    {
        mData->variant->applyCovInflation(*mData->upE, mData->covInflation, k);
    }

    mData->updateActive = true;
}



void EnsembleKalmanSmoother::assimilate(const Ref<const Array> z, 
                                        const ObservationOperator& H, const CovarianceOperator& R)
{
    ENDAS_ASSERT(mData->updateActive);

    // Nothing to do
    if (z.size() == 0) return;

    ENDAS_ASSERT(z.size() == H.nobs());

    auto& Eg = *mData->upE;
    int n = Eg.rows();
    int N = Eg.cols();

    // Global analysis
    if (mData->numDomains == 0)
    {
        // Compute whatever the variant will need from the global ensemble. This typically means 
        // some form of H(E)
        vector<Array2d> Egdata;
        {
            ENDAS_PERF_SCOPE(ProcessGlobalEnsemble);
            mData->variant->processGlobalEnsemble(Eg, H, mData->upK, Egdata);
        }

        // Compute the analysis ensemble transform
        {
            ENDAS_PERF_SCOPE(EnsembleTransform);
            if (!mData->upHaveX)
            {
                mData->variant->ensembleTransform(Eg, Egdata, z, R, mData->upK, mData->upX);
            }
            else
            {
                Matrix X(N, N);

                mData->variant->ensembleTransform(Eg, Egdata, z, R, mData->upK, X);
                mData->upX = mData->upX * X;
            }
        }

        // We only care about X if we are running smoother (if not, we will use upX conveniently 
        // as a pre-allocated array)
        mData->upHaveX = mData->upHaveX || mData->lag > 0;
    }
    else
    {
        ENDAS_NOT_IMPLEMENTED;
    }


}


//void 
//EnsembleKalmanSmoother::Data::assimilate(Ref<Array2d> E, void* Egdata, const Ref<const Array> z, 
//                                         const ObservationOperator& H, const CovarianceOperator& R)
//{
//
//}



void EnsembleKalmanSmoother::endAnalysis()
{
    ENDAS_ASSERT(mData->updateActive);

    // First the lagged smoother by updating all previous EnKF states with the transformation matrix
    // from the current analysis update. For lag=0 this does nothing

    mData->laggedSmoother(this->mOnResultFn, false);

    // Done with smoothing. Store the filter result for this update step and we're done

    // Filter only -> solution is available
    if (mData->lag == 0)
    {
        if (mOnResultFn) mOnResultFn(*mData->upE, mData->upK);
    }
    else
    {
        auto Ehandle = mData->cache->put(*mData->upE);
        mData->pushSmootherStep(mData->upK, Ehandle);
    }

    mData->updateActive = false;
}



void EnsembleKalmanSmoother::endSmoother()
{
    ENDAS_ASSERT(!mData->updateActive);
    if (mData->lag == 0) return;
    ENDAS_ASSERT(mOnResultFn);

    mData->laggedSmoother(this->mOnResultFn, true);
}




inline void smootherApplyForgetFactor(Ref<Matrix> X, double f)
{
    if (f == 1.0) return;

    // Apply forgetting factor `f` to `X`. This is done as X = ((X - I)*f ) + I
    for (int j = 0; j != X.cols(); j++)
    {
        for (int i = 0; i != X.rows(); i++)
        {
            if (i == j) X(i, j) = (X(i, j) - 1.0) * f + 1.0;
            else X(i, j)*= f;
        }
    }
    //X.diagonal().array() -= 1.0;  // X = X - I
    //X.array() *= forgettingFactor;
    //X.diagonal().array() += 1.0;  // X = X + I
}



void EnsembleKalmanSmoother::Data::laggedSmoother(OnResultFn& onresult, bool finishing)
{
    if (ksSteps.size() == 0) return;

    ENDAS_PERF_SCOPE(Smoother);

    int kend = ksSteps.size();
    for (int j = kend-1; j != kend-1-lag; j--)
    {
        if (j < 0) break;

        // Is the smoother result ready for this step?
        bool jIsResult = (j == kend - lag) || finishing;

        const auto& jdata = ksSteps[j];

        auto Ej = cache->get(jdata.E);
        ENDAS_ASSERT(Ej);


        // For global analysis, upE and upX are the global ensemble and transform arrays, so we only 
        // need to compute Ej * upX
        if (numDomains == 0)
        {
            if (upHaveX && !finishing)
            {
                smootherApplyForgetFactor(upX, forgettingFactor);
                Ej->array.matrix() = Ej->array.matrix() * upX;
            }
        }
        else
        {
            ENDAS_NOT_IMPLEMENTED;
        }
        
        // Have lagged result
        if (jIsResult)
        {
            if (onresult) onresult(Ej->array, jdata.k);

            // We will not need the ensemble for step k anymore
            cache->remove(jdata.E);
        }
        // We will still need Ej, let cache know it has been updated
        else
        {
            Ej->markDirty();
        }
    }

}





