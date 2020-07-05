#include <Endas/DA/SimpleObservationManager.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"


using namespace std;
using namespace endas;

struct SimpleObservationManager::Data
{
    Array obs; 
    Array2d obsCoords; 

    shared_ptr<const ObservationOperator> H;
    shared_ptr<const CovarianceOperator> R;

    const DomainPartitioning* partitioner;
    const TaperFn* taperFn;

    shared_ptr<const PartitionPointQuery> obsQuery;

    int currentDomain;
    int numDomains;

    Data(const Array o, Array2d oc, 
         shared_ptr<const ObservationOperator> h, shared_ptr<const CovarianceOperator> r)
    : obs(move(o)), obsCoords(move(oc)), H(h), R(r)
    {

    }
};



SimpleObservationManager::SimpleObservationManager(Array obs, Array2d obsCoords, 
                                                   shared_ptr<const ObservationOperator> H,
                                                   shared_ptr<const CovarianceOperator> R)
: mData(make_unique<Data>(move(obs), move(obsCoords), H, R))
{

}


SimpleObservationManager::~SimpleObservationManager()
{ }


void SimpleObservationManager::beginFetch(int k, const DomainPartitioning* partitioner, 
                                          const TaperFn* taperFn) const
{
    mData->partitioner = partitioner;
    mData->taperFn = taperFn;
    mData->currentDomain = 0;
    mData->numDomains = (partitioner)? partitioner->numLocalDomains() : 1;
    ENDAS_ASSERT(mData->numDomains >= 1);

    if (mData->numDomains > 1)
    {
        ENDAS_ASSERT(mData->obsCoords.size() > 0 || mData->obs.size() == 0);

        ENDAS_ASSERT(mData->taperFn);
        mData->obsQuery = partitioner->indexPoints(move(mData->obsCoords));
        ENDAS_ASSERT(mData->obsQuery);
    }
}


ObservationData SimpleObservationManager::fetchObservations() const
{
    // Global analysis
    if (mData->numDomains == 1)
    {
        return ObservationData(GlobalAnalysisDomainId, move(mData->obs), mData->H, mData->R); 
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

                return ObservationData(d, move(obsLocal), Hlocal, Rlocal);    
            }
        }

        return ObservationData();
    }

}
