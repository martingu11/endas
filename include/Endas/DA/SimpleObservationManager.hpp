/**
 * @file SimpleObservationManager.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_DA_SIMPLE_OBSERVATION_MANAGER_HPP__
#define __ENDAS_DA_SIMPLE_OBSERVATION_MANAGER_HPP__

#include <Endas/DA/ObservationManager.hpp>
#include <memory>

namespace endas
{

/** 
 * @addtogroup domain
 * @{ 
 */


/**
 * Simple observation manager.
 */
class ENDAS_DLL SimpleObservationManager  : public ObservationManager
{
public:

    SimpleObservationManager(const Ref<const Array> obs, const Ref<const Array2d> obsCoords, 
                             std::shared_ptr<const ObservationOperator> H,
                             std::shared_ptr<const CovarianceOperator> R);

    virtual ~SimpleObservationManager();

    virtual void beginFetch(int k, const DomainPartitioning* partitioner, 
                            const TaperFn* taperFn) const override;

    virtual ObservationData fetchObservations() const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};



/** @} */

}

#endif