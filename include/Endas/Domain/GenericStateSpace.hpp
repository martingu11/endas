/**
 * @file GenericStateSpace.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DOMAIN_GENERIC_STATESPACE_HPP__
#define __ENDAS_DOMAIN_GENERIC_STATESPACE_HPP__

#include <Endas/DA/StateSpace.hpp>
#include <Endas/DA/StateSpacePartitioning.hpp>
#include <Endas/Spatial/CoordinateSystem.hpp>

namespace endas
{

/** 
 * @addtogroup domain
 * @{ 
 */


/**
 * Trivial state space without assumptions on the structure of the state space.
 *
 * GenericStateSpace is likely only useful for low-dimensional toy problems where each state variable
 * can be updated independently. Both StateSpace and StateSpacePartitioning interfaces are 
 * implemented by this class as there is very little benefit in separating the two.
 * 
 * Under the implemented partitioning scheme, the index of each state variable is also its coordinate
 * and coordDim() is therefore 1. One-dimensional Euclidean distance is used by default to evaluate 
 * proximity although this can be replaced by user-defined function.
 */ 
class GenericStateSpace : public StateSpace, public StateSpacePartitioning
{
public:

    GenericStateSpace(index_t size);

    virtual index_t size() const override;
    virtual int coordDim() const override;    
    virtual const StateSpace& stateSpace() const override;

    virtual int numDomains() const override;
    virtual index_t getLocalStateSize(int d) const override;
    virtual void getLocalState(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const override;
    virtual void putLocalState(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const override;

    virtual std::shared_ptr<const PartitionPointQuery> 
    indexPoints(const Ref<const Array2d> coords) const override;

private:
    index_t mSize;
};



/** @} */

}

#endif