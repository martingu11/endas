/**
 * @file StateSpaceGeneric.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_STATESPACE_GENERIC_HPP__
#define __ENDAS_DA_STATESPACE_GENERIC_HPP__

#include "StateSpace.hpp"
#include <Endas/Spatial/CoordinateSystem.hpp>

#include <type_traits>


namespace endas
{


/**
 * Trivial state space with elements organized in one-dimensional sequence.
 * This is similar to a one-dimensional gridded state space but with simpler setup. 
 */ 
class Generic1dStateSpace : public StateSpace
{
public:

    Generic1dStateSpace(index_t size);

    virtual index_t size() const override;
    virtual const CoordinateSystem& crs() const override;

private:
    index_t mSize;
    EuclideanCS mCRS;
};


/**
 * State space partitioning scheme that assigns each state space variable into its own local domain.
 */ 
class GenericStateSpacePartitioning : public StateSpacePartitioning
{
public:

    /**
     * GenericStateSpacePartitioning constructor.
     * 
     * @param ss  State space. The instance is copied.
     */
    template <class SS, require_is_convertible<SS, const StateSpace> = true>
    GenericStateSpacePartitioning(const SS& ss)
    : GenericStateSpacePartitioning(std::make_shared<SS>(ss))
    { }

    /**
     * GenericStateSpacePartitioning constructor.
     * 
     * @param ss  State space.
     */
    GenericStateSpacePartitioning(std::shared_ptr<const StateSpace> ss); 

    virtual int numDomains() const override;
    virtual index_t getLocalStateSize(int d) const override;
    virtual void getLocalState(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const override;
    virtual void putLocalState(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const override;

private:
    
};



}

#endif