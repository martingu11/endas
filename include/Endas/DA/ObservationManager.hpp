/**
 * @file ObservationManager.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_DA_OBSERVATION_MANAGER_HPP__
#define __ENDAS_DA_OBSERVATION_MANAGER_HPP__

#include "ObservationOperator.hpp"
#include "CovarianceOperator.hpp"
#include "Domain.hpp"
#include "DomainPartitioning.hpp"
#include "Taper.hpp"

#include <memory>
#include <utility>


namespace endas
{

/** 
 * @addtogroup da
 * @{ 
 */


/** 
 * Special purpose local domain index specifying the global domain.
 * 
 * This should be the domain index returned from ObservationManager::fetchObservations()
 * when analysis localization is not in use.
 */
static constexpr int GlobalAnalysisDomainId = -1;


/**
 * Observation data returned by ObservationManager.
 * The data includes the observed values and the corresponding observation operator and 
 * observation error covariance operator.
 */
struct ObservationData
{
    /** Local domain index of the data or GlobalAnalysisDomainId. */
    int domain; 
    
    /** Array of observed values. */
    Array obs; 
    
    /** Observation operator corresponding to `obs`. */
    std::shared_ptr<const ObservationOperator> H;   
    
    /** Error covariance operator corresponding to `obs`. */
    std::shared_ptr<const CovarianceOperator> R;

    /** Constructs empty data entry. */
    ObservationData()
    : obs(emptyArray()), H(nullptr), R(nullptr)
    { }

    /** Constructs data entry. */
    ObservationData(int _domain, Array _obs, std::shared_ptr<const ObservationOperator> _H,
                    std::shared_ptr<const CovarianceOperator> _R)
    : domain(_domain), obs(std::move(_obs)), H(_H), R(_R)
    { }

    /** Returns `true` if the entry is empty. */
    bool empty() const { return obs.size() == 0; }
};




/**
 * Abstract observation manager.
 *
 * Observation managers are responsible for fetching observation data and the associated observation
 * operator and observation error covariance. 
 */
class ENDAS_DLL ObservationManager
{
public:

    virtual ~ObservationManager() { }

    /**
     * Called before assimilation of observations to notify ObservationManager that observations
     * will be requested. The manager should initialize its internal state so that observations 
     * can be fetched using fetchObservations().
     * 
     * @param k             Analysis time step index.
     * @param partitioner   State space partitioner if used (i.e. analysis is localized) or 
     *                      `nullptr`
     * @param taperFn       Observation covariance tapering function if used or `nullptr`
     */
    virtual void beginFetch(int k, const DomainPartitioning* partitioner, 
                            const TaperFn* taperFn) const = 0;



    /** 
     * Returns observation data.
     * 
     * If analysis is global (null partitioner was passed to beginFetch()), this should return all 
     * observations at once and set `domain` of the returned data instance to 
     * `GlobalAnalysisDomainId`. If analysis is local, this should return data for a single local 
     * domain and set `domain` to the domain index. The order in which data for domains is returned
     * can be arbitrary.
     * 
     * If there are no more observations, returns `Data()`.
     */
    virtual ObservationData fetchObservations() const = 0;

};


/** @} */

}

#endif