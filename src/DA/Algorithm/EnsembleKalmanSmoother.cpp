#include <Endas/DA/Algorithm/EnsembleKalmanSmoother.hpp>
#include <Endas/Core/Ensemble.hpp>
#include <Endas/Caching/MemoryArrayCache.hpp>
#include <Endas/Endas.hpp>

#include "../../Compatibility.hpp"

#include <vector>
#include <deque>
#include <iostream>


using namespace std;
using namespace endas;


using handle_t = ArrayCache::handle_t;

typedef SequentialEnsembleSmoother::OnResultFn OnResultFn;


//-------------------------------------------------------------------------------------------------
// EnKSVariant
//-------------------------------------------------------------------------------------------------

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


//-------------------------------------------------------------------------------------------------
// EnsembleKalmanSmoother
//-------------------------------------------------------------------------------------------------

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

    shared_ptr<ArrayCache> cache;

    // Data for the current update
    int upK;
    bool updateActive;
    bool upHaveAssimilatedObs;
   
    unique_ptr<Ref<Array2d>> upE; // Have to use pointer to Ref<> until we have Ptr<>
    Matrix upX; // Ensemble transform matrix
    bool upHaveX; // Is upX already holding a transform?

    // Smoother data
    double smForgetFactor;
    deque<ENKSData> ksSteps; // Data for all smoother steps

    // Data for localization
    shared_ptr<const DomainPartitioning> locSSP;
    shared_ptr<const TaperFn> locTaperFn;

    int locNumDomains;
    int locNumNonemptyDomains;
    index_t locTotalStateSize; // The sum of all locSSP.getLocalStateSize() values

    // Data for the augmented ensemble (global state including additional elements)
    Array2d locEaug;

    // (N*D) x N super-array of all X matrices  (D = number of nonempty domains)
    Array2d locX;

    // Array of flags whether locX is in use for a domain for ALL domains
    Eigen::Array<bool, Eigen::Dynamic, 1> locHaveX;
    
    // A 2xD array (here D is the number of ALL domains) containing the position where 
    // each domain starts in locEaug (row(0)) and the local state size (row(1))
    Eigen::Array<index_t, 2, Eigen::Dynamic> locStateLimits; 

    // Array containing the starting index of X matrices in locX of ALL domains
    Eigen::Array<index_t, 1, Eigen::Dynamic> locXLimits; 


    bool isLocalized() { return locNumDomains > 1; }


    Data() 
    : lag(0), covInflation(1.0), smForgetFactor(1.0), 
      updateActive(false), locNumDomains(0)
    { }

    Data(const Data&) = delete;

    void pushSmootherStep(int k, handle_t E)
    {
        ksSteps.push_back({ k, E });
    }

    template <class Fn>
    void foreachNonemptyDomain(Fn fn)
    {
        int i = 0;
        for (int d = 0; d != locNumDomains; d++)
        {
            index_t start = locStateLimits.col(d)(0);
            index_t nloc = locStateLimits.col(d)(1);
            if (nloc == 0) continue;
            fn(d, i++, start, nloc);
        }
    }

    void assimilateOne(const ObservationData& odata, 
                       const Ref<const Array2d> Eg, Ref<Array2d> E, Ref<Matrix> X, bool& haveX);

    void laggedSmoother(OnResultFn& onresult, bool finishing);

    
    // Unpacks E into the augmented ensemble (locEaug) when local domains are in use.
    void unpackEnsemble(const Ref<const Array2d> E, Ref<Array2d> Eaug);

    // Packs the augmented ensemble (locEaug) back to E when local domains are in use.
    void packEnsemble(const Ref<const Array2d> Eaug, Ref<Array2d> E);

};



EnsembleKalmanSmoother::EnsembleKalmanSmoother(const EnKSVariant& variant, int n, int N,
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
    mData->smForgetFactor = factor;
}


void EnsembleKalmanSmoother::localize(shared_ptr<const DomainPartitioning> partitioner, 
                                      shared_ptr<const TaperFn> taperFn)
{
    ENDAS_ASSERT(partitioner);
    mData->locSSP = partitioner;
    mData->locTaperFn = taperFn;

    mData->locNumDomains = mData->locSSP->numLocalDomains();
    ENDAS_ASSERT(mData->locNumDomains > 0);

    // Effectively global analysis
    if (mData->locNumDomains == 1)
    {
        globalize();
        return;
    }

    // Go once through all domains and check their state size. It is possible that the overall state 
    // size (locTotalStateSize) is larger than n if local domains contains some padding. Therefore we 
    // will use an 'augmented' global state with additional elements, if needed. 

    mData->locNumNonemptyDomains = 0;
    mData->locTotalStateSize = 0;
    mData->locStateLimits.resize(2, mData->locNumDomains);
    mData->locXLimits.resize(1, mData->locNumDomains);

    for (int d = 0; d != mData->locNumDomains; d++)
    {
        index_t nloc = mData->locSSP->getLocalSize(d);
        ENDAS_ASSERT(nloc >= 0); // 0 is allowed (empty domain) as it simplifies SSP implementations

        mData->locStateLimits.col(d)(0) = mData->locTotalStateSize;
        mData->locStateLimits.col(d)(1) = nloc;
        mData->locTotalStateSize+= nloc;

        mData->locXLimits(d) = mData->locNumNonemptyDomains * mData->N;

        if (nloc > 0) ++mData->locNumNonemptyDomains;
    }

    // Pre-allocate the augmented global ensemble and X arrays. 
    mData->locEaug.resize(mData->locTotalStateSize, mData->N);
    mData->locX.resize(mData->locNumNonemptyDomains * mData->N, mData->N);
    mData->locHaveX.resize(mData->locNumDomains * mData->N, 1);
}


void EnsembleKalmanSmoother::globalize()
{
    mData->locSSP.reset();
    mData->locTaperFn.reset();
    mData->locNumDomains = 1; 
    mData->locTotalStateSize = mData->n;
    mData->locStateLimits.resize(2, 0);
    mData->locXLimits.resize(1, 0);
    mData->locEaug.resize(0, 0);
    mData->locX.resize(0, 0);  
    mData->locHaveX.resize(0, 0);
}


void EnsembleKalmanSmoother::Data::unpackEnsemble(const Ref<const Array2d> E, Ref<Array2d> Eaug)
{
    foreachNonemptyDomain([&](int d, int i, index_t start, index_t nloc)
    {
        locSSP->getLocal(d, E, Eaug.block(start, 0, nloc, N));
    });
    /*for (int d = 0; d != locNumDomains; d++)
    {
        index_t start = locStateLimits.col(d)(0);
        index_t nloc = locStateLimits.col(d)(1);
        if (nloc == 0) continue;
        
    }*/
}

void EnsembleKalmanSmoother::Data::packEnsemble(const Ref<const Array2d> Eaug, Ref<Array2d> E)
{
    foreachNonemptyDomain([&](int d, int i, index_t start, index_t nloc)
    {
        locSSP->putLocal(d, Eaug.block(start, 0, nloc, N), E);
    });

    /*for (int d = 0; d != locNumDomains; d++)
    {
        index_t start = locStateLimits.col(d)(0);
        index_t nloc = locStateLimits.col(d)(1);
        if (nloc == 0) continue;
        locSSP->putLocalState(d, Eaug.block(start, 0, nloc, N), E);
    }*/
}


void EnsembleKalmanSmoother::beginSmoother(const Ref<const Array2d> E0, int k0)
{
    ENDAS_ASSERT(!mData->updateActive);
    ENDAS_ASSERT(E0.rows() == mData->n);
    ENDAS_ASSERT(E0.cols() == mData->N);

    if (mData->lag == 0) return;

    mData->ksSteps.clear();
    mData->upE.reset();

    handle_t E0handle;
    if (!mData->isLocalized())
    {
         E0handle = mData->cache->put(E0);
    }
    else
    {
        mData->unpackEnsemble(E0, mData->locEaug);
        E0handle = mData->cache->put(mData->locEaug);
    }
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
    mData->upHaveAssimilatedObs = false;

    // Apply covariance inflation 
    if (mData->covInflation != 1.0)
    {
        mData->variant->applyCovInflation(*mData->upE, mData->covInflation, k);
    }

    if (mData->isLocalized())
    {
        mData->locHaveX.fill(false);

        mData->unpackEnsemble(E, mData->locEaug);
    }

    mData->updateActive = true;
}



void EnsembleKalmanSmoother::assimilate(const ObservationManager& omgr)
{
    ENDAS_ASSERT(mData->updateActive);

    ENDAS_PERF_SCOPE(AssimilateObservations);

    auto& Eg = *mData->upE;
    int N = Eg.cols();

    omgr.beginFetch(mData->upK, mData->locSSP.get(), mData->locTaperFn.get());


    // Global analysis
    if (!mData->isLocalized())
    {
        ObservationData odata = omgr.fetchObservations();
        if (!odata.empty())
        {
            ENDAS_ASSERT(odata.domain == GlobalAnalysisDomainId);
            mData->assimilateOne(odata, Eg, Eg, mData->upX, mData->upHaveX);
        }
    }
    // Localized analysis
    else
    {
        // Have assimilated observations already in this analysis step -> need to reconstruct
        // global ensemble to be able to call H() on it
        if (mData->upHaveAssimilatedObs)
        {
            mData->packEnsemble(mData->locEaug, Eg);
        }

        //cout << "Local assim " << mData->locNumDomains << endl;

        // Fetch observations from the manager and assimilate...
        ObservationData odata;
        while ((odata = omgr.fetchObservations()).empty() == false)
        {
            int d = odata.domain;

            //cout << "Got data for " << d << endl;

            ENDAS_ASSERT(d != GlobalAnalysisDomainId && d >= 0 && d < mData->locNumDomains);

            index_t start = mData->locStateLimits.col(d)(0);
            index_t nloc = mData->locStateLimits.col(d)(1);
            index_t Xstart = mData->locXLimits(d);

            // Local (augmented) state and X
            auto Ed = mData->locEaug.block(start, 0, nloc, N);
            auto Xd  = mData->locX.block(Xstart, 0, N, N);
            bool& haveXd = mData->locHaveX(d);

            //cout << "  Assim " << d << endl;

            mData->assimilateOne(odata, Eg, Ed, Xd.matrix(), haveXd);

            //cout << "  Done " << d << endl;
        }
    }
}


void EnsembleKalmanSmoother::Data::assimilateOne(const ObservationData& odata, 
                                                 const Ref<const Array2d> Eg, Ref<Array2d> E, 
                                                 Ref<Matrix> X, bool& haveX)
{
    ENDAS_ASSERT(odata.H);
    ENDAS_ASSERT(odata.R);
    ENDAS_ASSERT(odata.H->nobs() == odata.obs.size());
    ENDAS_ASSERT(odata.R->size() == odata.obs.size());


    // Compute whatever the variant will need from the global ensemble. This typically means 
    // some form of H(E)
    vector<Array2d> Egdata;
    {
        ENDAS_PERF_SCOPE(ProcessGlobalEnsemble);
        variant->processGlobalEnsemble(Eg, *odata.H, upK, Egdata);
    }

    // Compute the analysis ensemble transform
    {
        ENDAS_PERF_SCOPE(EnsembleTransform);
        if (!haveX)
        {
            variant->ensembleTransform(E, Egdata, odata.obs, *odata.R, upK, X);
        }
        else
        {
            Matrix XX(N, N);
            variant->ensembleTransform(E, Egdata, odata.obs, *odata.R, upK, XX);
            X = X * XX;
        }
    }

    // We only care about X if we are running smoother (if not, we will use upX conveniently 
    // as a pre-allocated array)
    haveX = haveX || lag > 0;  

    upHaveAssimilatedObs = true; 
}


void EnsembleKalmanSmoother::endAnalysis()
{
    ENDAS_ASSERT(mData->updateActive);

    // First the lagged smoother by updating all previous EnKF states with the transformation matrix
    // from the current analysis update. For lag=0 this does nothing

    mData->laggedSmoother(this->mOnResultFn, false);

    // Done with smoothing. Store the filter result for this update step and we're done. 
    // If localized analysis was done and `upE` is not up to date, we need to reconstruct it first
    if (mData->isLocalized() && (mData->upHaveAssimilatedObs || mData->lag > 0))
    {
        mData->packEnsemble(mData->locEaug, *mData->upE);
    }


    // Filter only -> solution is available
    if (mData->lag == 0)
    {
        if (mOnResultFn) mOnResultFn(*mData->upE, mData->upK);
    }
    else
    {
        handle_t Ehandle = (!mData->isLocalized())? 
            mData->cache->put(*mData->upE) :
            mData->cache->put(mData->locEaug);

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

        if (!finishing)
        {
            // For global analysis, upE and upX are the global ensemble and transform arrays, so we only 
            // need to compute Ej * upX
            if (!isLocalized())
            {
                if (upHaveX)
                {
                    smootherApplyForgetFactor(upX, smForgetFactor);
                    Ej->array.matrix() = Ej->array.matrix() * upX;
                }
            }
            // For local analyses we operate on the augmented ensemble and use the per-domain X matrices
            else
            {
                foreachNonemptyDomain([&](int d, int i, index_t start, index_t nloc)
                {
                    if (!locHaveX(d)) return;

                    // Local (augmented) state and X
                    auto Ed = Ej->array.block(start, 0, nloc, N);
                    auto Xd  = locX.block(i*N, 0, N, N).matrix();
                    smootherApplyForgetFactor(Xd, smForgetFactor);

                    Ed.matrix() = Ed.matrix() * Xd;
                });
            }
        }
        
        // Have lagged result
        if (jIsResult)
        {
            if (onresult) 
            {
                // Localized analysis -> result is in Ej which is the augmented ensemble. Reconstruct
                // original first. Here we resue `upE` since it only holds useful data if no observations
                // were assimilated at this step. In that case we will have to reconstruct `upE` from the 
                // augmented ensemble later (in endAnalysis()) at some cost.
                if (isLocalized())
                {   
                    packEnsemble(Ej->array, *upE);
                    onresult(*upE, jdata.k);
                }
                else 
                    onresult(Ej->array, jdata.k);
            }

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

