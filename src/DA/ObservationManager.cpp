
#include <Endas/DA/ObservationManager.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;



ObservationManager::~ObservationManager()
{ }



SimpleObservationManager::SimpleObservationManager(const Ref<const Array> obs, 
                                                   shared_ptr<const ObservationOperator> H,
                                                   shared_ptr<const CovarianceOperator> R)
: mObs(obs), mH(H), mR(R)
{ }


ObservationManager::Data SimpleObservationManager::getObservations(int domain) const
{
    if (domain == NullDomain)
    {
        return Data(mObs, mH, mR);
    }
    else
    {
        return Data();
    }

}