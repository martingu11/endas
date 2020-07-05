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
 * 
 * This is a very basic implementation of the ObservationManager protocol that simply serves
 * the given observation array, coordinates, observation operator and error covariance. 
 * The manager supports localized analysis and relies on point indexing capabilities of the 
 * domain partitioning scheme as well as ObservationOperator::subset() and 
 * CovarianceOperator::subset(). 
 * 
 * The implementation may be a good start for learning how to implement more complex observation 
 * managers. An example may be an observation manager that fetches observations from a web data
 * interface or a database rather than relying on all observations to be available up-front.
 * 
 * .. attention::
 *    The manager can only be used to provide observations for a single data assimilation step
 *    assimilating the same observations multiple times makes very little sense anyway.
 *    Therefore, the intended usage is to construct new manager each time observations are 
 *    assimilated. This allows the observation and coordinate arrays top be moved rather than
 *    copied during the beginFetch() and fetchObservations() calls. 
 */
class ENDAS_DLL SimpleObservationManager  : public ObservationManager
{
public:

    SimpleObservationManager(Array obs, Array2d obsCoords, 
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