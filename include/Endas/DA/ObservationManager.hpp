/**
 * @file ObservationManager.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_DA_OBSERVATION_MANAGER_HPP__
#define __ENDAS_DA_OBSERVATION_MANAGER_HPP__

#include "ObservationOperator.hpp"
#include "CovarianceOperator.hpp"
#include "StateSpace.hpp"
#include "Taper.hpp"

#include <memory>


namespace endas
{

/** 
 * Special purpose local domain index specifying "no domain".
 */
static constexpr int NullDomain = -1;


/**
 * Abstract observation manager.
 *
 * Observation managers are responsible for fetching observation data and the associated observation
 * operator and observation error covariance. 
 */
class ENDAS_DLL ObservationManager
{
public:

    /**
     * Observation data returned by ObservationManager.
     * The data includes the observed values and the corresponding observation operator and 
     * observation error covariance operator.
     */
    struct Data
    {
        Ref<const Array> obs;
        std::shared_ptr<const ObservationOperator> H;
        std::shared_ptr<const CovarianceOperator> R;

        Data()
        : obs(emptyArray()), H(nullptr), R(nullptr)
        { }

        Data(Ref<const Array> _obs, std::shared_ptr<const ObservationOperator> _H,
             std::shared_ptr<const CovarianceOperator> _R)
        : obs(_obs), H(_H), R(_R)
        { }
    };

    virtual ~ObservationManager();


    //virtual void beginAnalysis(int k, std::shared_ptr<const StateSpacePartitioning> partitioner,
    //                           const TaperFn* taperFn) const = 0;



    /** 
     * Returns observational data for given analysis domain.
     * 
     * @param domain       Index of the analysis domain for which data is to be retrieved. 
     *                     endas::NullDomain is passed if analysis is global.
     * @param partitioner  
     */
    virtual Data getObservations(int domain) const = 0;

};



class ENDAS_DLL SimpleObservationManager  : public ObservationManager
{
public:

    SimpleObservationManager(const Ref<const Array> obs, 
                             std::shared_ptr<const ObservationOperator> H,
                             std::shared_ptr<const CovarianceOperator> R);

                      

    virtual Data getObservations(int domain) const override;

private:
    
    const Ref<const Array> mObs; 
    std::shared_ptr<const ObservationOperator> mH;
    std::shared_ptr<const CovarianceOperator> mR;
};





}

#endif