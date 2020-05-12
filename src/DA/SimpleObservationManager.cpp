#include <Endas/DA/SimpleObservationManager.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"


using namespace std;
using namespace endas;

struct SimpleObservationManager::Data
{
    const Ref<const Array> obs; 
    const Ref<const Array2d> obsCoords; 

    shared_ptr<const ObservationOperator> H;
    shared_ptr<const CovarianceOperator> R;

    const StateSpacePartitioning* partitioner;
    const TaperFn* taperFn;

    shared_ptr<const PartitionPointQuery> obsQuery;

    int currentDomain;
    int numDomains;

    Data(const Ref<const Array> o, const Ref<const Array2d> oc, 
         shared_ptr<const ObservationOperator> h, shared_ptr<const CovarianceOperator> r)
    : obs(o), obsCoords(oc), H(h), R(r)
    {

    }
};



SimpleObservationManager::SimpleObservationManager(const Ref<const Array> obs, 
                                                   const Ref<const Array2d> obsCoords, 
                                                   shared_ptr<const ObservationOperator> H,
                                                   shared_ptr<const CovarianceOperator> R)
: mData(make_unique<Data>(obs, obsCoords, H, R))
{

}


SimpleObservationManager::~SimpleObservationManager()
{ }


void SimpleObservationManager::beginFetch(int k, const StateSpacePartitioning* partitioner, 
                                          const TaperFn* taperFn) const
{
    mData->partitioner = partitioner;
    mData->taperFn = taperFn;
    mData->currentDomain = 0;
    mData->numDomains = (partitioner)? partitioner->numDomains() : 1;
    ENDAS_ASSERT(mData->numDomains >= 1);

    if (mData->numDomains > 1)
    {
        ENDAS_ASSERT(mData->taperFn);
        mData->obsQuery = partitioner->indexPoints(mData->obsCoords);
        ENDAS_ASSERT(mData->obsQuery);
    }
}


ObservationManager::Data SimpleObservationManager::fetchObservations() const
{

    // Global analysis
    if (mData->numDomains == 1)
    {
        /// @todo This implies a copy of obs!!!
        return ObservationManager::Data(GlobalAnalysisDomainId, mData->obs, mData->H, mData->R); 
    }
    // Local analysis
    else
    {
        IndexArray obsIndices;

        // Cycle through domains until we have one with observations
        while (mData->currentDomain != mData->numDomains)
        {
            // Query which observations we will need for this domain. This will be all that fall within the
            // localization taper function support range.
            obsIndices.clear();
            mData->obsQuery->rangeQuery(mData->currentDomain, mData->taperFn->supportRange(), obsIndices);

            int d = mData->currentDomain++;

            if (obsIndices.size() > 0)
            {
                Array obsLocal(obsIndices.size());
                select(mData->obs, obsIndices, obsLocal);

                auto Hlocal = mData->H->subset(obsIndices);
                auto Rlocal = mData->R->subset(obsIndices);
                ENDAS_ASSERT(Hlocal);
                ENDAS_ASSERT(Rlocal);

                return ObservationManager::Data(d, move(obsLocal), Hlocal, Rlocal);    
            }
        }

        return ObservationManager::Data();
    }

}
