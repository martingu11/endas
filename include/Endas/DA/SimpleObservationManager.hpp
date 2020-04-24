/**
 * @file SimpleObservationManager.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_DA_SIMPLE_OBSERVATION_MANAGER_HPP__
#define __ENDAS_DA_SIMPLE_OBSERVATION_MANAGER_HPP__

#include "ObservationManager.hpp"
#include <memory>

namespace endas
{



/**
 * Basic observation manager.
 */
class ENDAS_DLL SimpleObservationManager  : public ObservationManager
{
public:

    SimpleObservationManager(const Ref<const Array> obs, const Ref<const Array2d> obsCoords, 
                             std::shared_ptr<const ObservationOperator> H,
                             std::shared_ptr<const CovarianceOperator> R);

    virtual ~SimpleObservationManager();

    virtual void beginFetch(int k, const StateSpacePartitioning* partitioner, 
                            const TaperFn* taperFn) const override;

    virtual Data fetchObservations() const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};





}

#endif