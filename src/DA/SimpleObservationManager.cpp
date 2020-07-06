#include <Endas/DA/SimpleObservationManager.hpp>
#include <Endas/DA/Taper.hpp>
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
        // We will do covariance tapering if the taper function is not NoTaper and has 
        // support range > 0
        bool doCovTapering = 
            dynamic_cast<const NoTaper*>(mData->taperFn) == nullptr &&
            mData->taperFn->supportRange() > 0;

        IndexArray obsIndices;
        PartitionPointQuery::DistanceArray obsDistances;

        // Cycle through domains until we have one with observations
        while (mData->currentDomain != mData->numDomains)
        {
            // Query which observations we will need for this domain. This will be all that fall within the
            // localization taper function support range.
            obsIndices.clear();
            obsDistances.clear();

            mData->obsQuery->rangeQuery(
                mData->currentDomain, mData->taperFn->supportRange(), 
                obsIndices, 
                (doCovTapering)? &obsDistances : nullptr);

            if (doCovTapering)
            {
                ENDAS_ASSERT(obsDistances.size() == obsIndices.size());
            }

            int d = mData->currentDomain++;

            if (obsIndices.size() > 0)
            {
                Array obsLocal(obsIndices.size());
                select(mData->obs, obsIndices, obsLocal);

                auto Hlocal = mData->H->subset(obsIndices);
                auto Rlocal = mData->R->subset(obsIndices);
                ENDAS_ASSERT(Hlocal);
                ENDAS_ASSERT(Rlocal);
                ENDAS_ASSERT(Hlocal->nobs() == obsIndices.size());
                ENDAS_ASSERT(Rlocal->size() == obsIndices.size());

                // Taper the observation error covariance matrix by the given tapering function.
                // In practice, this means multiplying the covariance entries so that (co)varinace of 
                // observations far away is increaed, effectively limiting their influence. 
                // Currently we only handle diagonal covariance matrices here and set
                //   
                //    cov(i,i) = cov(i,i) * 1 / taper(distance(i, domain))
                // 
                // where cov(i,i) is the diagonal covariance element, distance(i, domain) is the 
                // distance of i-th observation to the domain and taper(d) is the taper function 
                // value dor distance d. 
                
                // Todo: Remove observations for which the tapering factor is very close to zero to avoid 
                // numerical issues with the covariance!

                const DiagonalCovariance* RasDiag = dynamic_cast<const DiagonalCovariance*>(Rlocal.get());
                if (doCovTapering && RasDiag)
                {
                    // A copy of the inverse diagonal. The formula above is implemented by tapering the 
                    // inverse diagonal to avoid numerical issues. Most filters work with inverse of the
                    // covariance anyway.

                    Array RdiagInv = RasDiag->inverseDiagonal(); 
                    ENDAS_ASSERT(RdiagInv.size() == obsDistances.size());

                    // Taper the inverse diagonal entries by taper(distance(i, domain))
                    Eigen::Map<Array> distancesAsArray(obsDistances.data(), obsDistances.size());
                    mData->taperFn->taper(RdiagInv, distancesAsArray, RdiagInv);

                    // Construct new covariance operator using the tapered inverse diagonal
                    Rlocal = make_shared<DiagonalCovariance>(move(RdiagInv), true);

                }

                return ObservationData(d, move(obsLocal), Hlocal, Rlocal);    
            }
        }

        return ObservationData();
    }

}
