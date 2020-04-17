#include <Endas/Algorithm/EnsembleKalmanSmoother.hpp>
#include <Endas/Caching/MemoryArrayCache.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"

#include <vector>


using namespace std;
using namespace endas;

using handle_t = ArrayCache::handle_t;


/* Smoother data for one time step. */
struct ENKSData
{
    int k;
    handle_t xf, Pf, xa, Pa;
};



struct EnsembleKalmanSmoother::Data
{
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
    vector<ENKSData> ksSteps; // Data for all smoother steps

    Data() : lag(0), updateActive(false) { }

    Data(const Data&) = delete;


    void pushSmootherStep(int k, handle_t xf, handle_t Pf, handle_t xa, handle_t Pa)
    {
        ksSteps.push_back({ k, xf, Pf, xa, Pa });
    }

    ENKSData popSmootherStep()
    {
        ENDAS_ASSERT(ksSteps.size() > 0);
        auto ksdata = ksSteps.back();
        ksSteps.pop_back();
        return ksdata;
    }


};


EnsembleKalmanSmoother::EnsembleKalmanSmoother(const EnKFVariant& variant, int lag, 
                                               shared_ptr<ArrayCache> cache)
: mData()
{
    ENDAS_ASSERT(lag >= 0);

    mData->lag = lag;
    mData->cache = (cache)? cache : make_shared<MemoryArrayCache>();
}


EnsembleKalmanSmoother::~EnsembleKalmanSmoother()
{ }




void EnsembleKalmanSmoother::beginSmoother(const Ref<const Array2d>E0, int k0)
{
    ENDAS_ASSERT(!mData->updateActive);
    if (mData->lag == 0) return;

    mData->ksSteps.clear();
}



void EnsembleKalmanSmoother::beginAnalysis(Ref<Array2d> E, int k)
{
    ENDAS_ASSERT(!mData->updateActive);

    /*mData->upX = make_unique<Ref<ColVec>>(x);
    mData->upP = make_unique<Ref<Matrix>>(P);
    mData->upK = k;
    mData->smXf = mData->cache->put(x);
    mData->smPf = mData->cache->put(P);*/

    mData->updateActive = true;
}

void EnsembleKalmanSmoother::assimilate(const Ref<const Array> z, const ObservationOperator& H, 
                                        const CovarianceOperator& R)
{
    ENDAS_ASSERT(mData->updateActive);

    // Nothing to do
    if (z.size() == 0) return;

    /*ENDAS_ASSERT(z.size() == H.rows());

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
    upP.noalias() -= PHt * cholF.solve(H * upP);*/
}


void EnsembleKalmanSmoother::endAnalysis()
{
    ENDAS_ASSERT(mData->updateActive);

    // Filter only -> solution is available
    /*if (mData->lag == 0)
    {
        if (mOnResultFn) mOnResultFn(*mData->upX, *mData->upP, mData->upK);
    }
    // Smoother -> store smoother data for this time step
    else
    {
        auto xhandle = mData->cache->put(*mData->upX);
        auto Phandle = mData->cache->put(*mData->upP);
        mData->pushSmootherStep(mData->upK, mData->smXf, mData->smPf, xhandle, Phandle);
    }*/

    mData->updateActive = false;
}



void EnsembleKalmanSmoother::endSmoother()
{
    ENDAS_ASSERT(!mData->updateActive);
    if (mData->lag == 0) return;
    if (mData->ksSteps.size() == 0) return;

    ENDAS_ASSERT(mOnResultFn);

    /*
    // The last filter solution is also the smoother solution
    auto sm_last = mData->popSmootherStep();

    // We won't need model tangent linear for the last time step
    mData->model.stepFinished(sm_last.k);

    auto xs = mData->cache->pop(sm_last.xa);
    auto Ps = mData->cache->pop(sm_last.Pa);
    ENDAS_ASSERT(xs);
    ENDAS_ASSERT(Ps);

    mData->onResultFn(*xs, *Ps, sm_last.k);

    // Backward recursion 
    Matrix MtPa(Ps->rows(), Ps->cols());
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

        MtPa = *Pa;
        mData->model.tl(MtPa, sm_k.k);

        Eigen::LLT<Matrix> cholPf(*Pf);
        Matrix J = cholPf.solve(MtPa).transpose();
        
        //if (mData->forgetFactor != 1.0) J*= mData->forgetFactor;

        xa->matrix().noalias() += J * (xs->matrix() - xf->matrix());
        Pa->matrix().noalias() += J * (Ps->matrix() - Pf->matrix()).transpose() * J.transpose();
        //Pa->matrix().noalias() += J * (J * (Ps->matrix() - Pf->matrix()).transpose());

        xs = xa;
        Ps = Pa;
        sm_k1 = sm_k;

        mData->onResultFn(*xs, *Ps, sm_k.k);
    }*/
}




