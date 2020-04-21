#include <Endas/Algorithm/KalmanSmoother.hpp>
#include <Endas/Caching/MemoryArrayCache.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"

#include <vector>


using namespace std;
using namespace endas;

using handle_t = ArrayCache::handle_t;


/* Smoother data for one time step. */
struct KSData
{
    int k;
    handle_t xf, Pf, xa, Pa;
};



struct KalmanSmoother::Data
{
    const Model& model;
    int lag;

    shared_ptr<ArrayCache> cache;

    // Data for the current update
    bool updateActive;
    unique_ptr<Ref<ColVec>> upX; // Have to use pointer to Ref<> until we have Ptr<>
    unique_ptr<Ref<Matrix>> upP;
    int upK;

    // Smoother data
    handle_t smXf;  // Handle to the prior state vector for current update
    handle_t smPf;  // Handle to the prior error cov. for current update
    vector<KSData> ksSteps; // Data for all smoother steps

    Data(const Model& m)
    : model(m), lag(0),
      updateActive(false)
    { }

    Data(const Data&) = delete;


    void pushSmootherStep(int k, handle_t xf, handle_t Pf, handle_t xa, handle_t Pa)
    {
        ksSteps.push_back({ k, xf, Pf, xa, Pa });
    }

    KSData popSmootherStep()
    {
        ENDAS_ASSERT(ksSteps.size() > 0);
        auto ksdata = ksSteps.back();
        ksSteps.pop_back();
        return ksdata;
    }


};


KalmanSmoother::KalmanSmoother(const Model& model, int lag, shared_ptr<ArrayCache> cache)
: mData(make_unique<Data>(model))
{
    ENDAS_ASSERT(lag >= 0);

    mData->lag = lag;
    mData->cache = (cache)? cache : make_shared<MemoryArrayCache>();
}


KalmanSmoother::~KalmanSmoother()
{ }


void KalmanSmoother::beginSmoother(const Ref<const Array>x0, const Ref<const Matrix> P0, int k0)
{
    ENDAS_ASSERT(!mData->updateActive);
    if (mData->lag == 0) return;

    mData->ksSteps.clear();
    mData->upX.reset();
    mData->upP.reset();

    auto x0handle = mData->cache->put(x0);
    auto P0handle = mData->cache->put(P0);
    mData->pushSmootherStep(k0, ArrayCache::NullHandle, ArrayCache::NullHandle, x0handle, P0handle);
}


void KalmanSmoother::forecast(Ref<Array> x, Ref<Matrix> P, const Ref<const Matrix> Q, int k, double dt)
{
    ENDAS_ASSERT(!mData->updateActive);
    ENDAS_PERF_SCOPE(Forecast);

    // Propagate model 
    {
        ENDAS_PERF_SCOPE(Model);
        mData->model.apply(x, k, dt, true);
    }

    // Propagate error covariance Pk+1 = M Pk M' + Q
    {
        ENDAS_PERF_SCOPE(ModelTangentLinear);
        mData->model.tl(P, k);
    }
    {
        ENDAS_PERF_SCOPE(ModelAdjoint);
        mData->model.adj(P, k);
    }

    if (Q.size() > 0) P+= Q;

    // Model can drop any data for step `k`
    if (mData->lag == 0) mData->model.stepFinished(k);
}


void KalmanSmoother::beginAnalysis(Ref<Array> x, Ref<Matrix> P, int k)
{
    ENDAS_ASSERT(!mData->updateActive);

    mData->upX = make_unique<Ref<ColVec>>(x);
    mData->upP = make_unique<Ref<Matrix>>(P);
    mData->upK = k;
    mData->smXf = mData->cache->put(x);
    mData->smPf = mData->cache->put(P);

    mData->updateActive = true;
}

void KalmanSmoother::assimilate(const Ref<const Array> z, const Ref<const Matrix> H, const Ref<const Matrix> R)
{
    ENDAS_ASSERT(mData->updateActive);

    ENDAS_PERF_SCOPE(Update);

    // Nothing to do
    if (z.size() == 0) return;

    ENDAS_ASSERT(z.size() == H.rows());

    auto& upX = *mData->upX;
    auto& upP = *mData->upP;

    Matrix PHt = upP * H.transpose();

    // Innovation    
    Matrix dz = z;
    dz.noalias() -= H * upX;

    // Innovation covariance F = H P H' + R
    Matrix F = R;
    F.noalias() += H * PHt;

    // State update as x = x + P H' F^-1 dz
    Eigen::LLT<Ref<Matrix>> cholF(F);
    upX.noalias() +=  PHt * cholF.solve(dz);

    // Error covariance update as P = P - P H' F^-1 H P 
    upP.noalias() -= PHt * cholF.solve(H * upP);
}


void KalmanSmoother::endAnalysis()
{
    ENDAS_ASSERT(mData->updateActive);

    // Filter only -> solution is available
    if (mData->lag == 0)
    {
        if (mOnResultFn) mOnResultFn(*mData->upX, *mData->upP, mData->upK);
    }
    // Smoother -> store smoother data for this time step
    else
    {
        auto xhandle = mData->cache->put(*mData->upX);
        auto Phandle = mData->cache->put(*mData->upP);
        mData->pushSmootherStep(mData->upK, mData->smXf, mData->smPf, xhandle, Phandle);
    }

    mData->updateActive = false;
}



void KalmanSmoother::endSmoother()
{
    ENDAS_ASSERT(!mData->updateActive);
    if (mData->lag == 0) return;
    if (mData->ksSteps.size() == 0) return;

    ENDAS_ASSERT(mOnResultFn);

    ENDAS_PERF_SCOPE(Smoother);


    /// @todo Currently we're having RTS implementation only but fixed-lag KS would 
    /// be desirable too.

    // The last filter solution is also the smoother solution
    auto sm_last = mData->popSmootherStep();

    // We won't need model tangent linear for the last time step
    mData->model.stepFinished(sm_last.k);

    auto xs = mData->cache->pop(sm_last.xa);
    auto Ps = mData->cache->pop(sm_last.Pa);
    ENDAS_ASSERT(xs);
    ENDAS_ASSERT(Ps);

    mOnResultFn(xs->array, Ps->array, sm_last.k);

    // Backward recursion 
    Matrix MtPa(Ps->array.rows(), Ps->array.cols());
    auto sm_k1 = sm_last;

    while (mData->ksSteps.size() > 0)
    {
        auto sm_k = mData->popSmootherStep();

        auto xf = mData->cache->pop(sm_k1.xf);
        auto Pf = mData->cache->pop(sm_k1.Pf);
        auto xa = mData->cache->pop(sm_k.xa);
        auto Pa = mData->cache->pop(sm_k.Pa);
        ENDAS_ASSERT(xf);
        ENDAS_ASSERT(Pf);
        ENDAS_ASSERT(xa);
        ENDAS_ASSERT(Pa);

        MtPa = Pa->array;
        mData->model.tl(MtPa, sm_k.k);

        Eigen::LLT<Matrix> cholPf(Pf->array);
        Matrix J = cholPf.solve(MtPa).transpose();
        
        //if (mData->forgetFactor != 1.0) J*= mData->forgetFactor;

        xa->array.matrix().noalias() += J * (xs->array.matrix() - xf->array.matrix());
        Pa->array.matrix().noalias() += J * (Ps->array.matrix() - Pf->array.matrix()).transpose() * J.transpose();
        //Pa->matrix().noalias() += J * (J * (Ps->matrix() - Pf->matrix()).transpose());

        xs = xa;
        Ps = Pa;
        sm_k1 = sm_k;

        mOnResultFn(xs->array, Ps->array, sm_k.k);
    }
}




