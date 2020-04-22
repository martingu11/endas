/**
 * @file StateSpace.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_STATESPACE_HPP__
#define __ENDAS_DA_STATESPACE_HPP__

#include <Endas/Core/Core.hpp>
#include <Endas/Spatial/Spatial.hpp>
#include <Endas/Spatial/CoordinateSystem.hpp>


namespace endas
{


/**
 * Base class for state space implementations.
 * 
 * This is a generic state space abstraction that only provides information about the size of the
 * state space and the coordinate system for measuring proximity within the space.
 * 
 * Please note that StateSpace or derived classes are only needed whe using analysis localization 
 * or with observation managers that spatially map observations to state variables. In other cases 
 * no information needs to be known about the state space.
 */ 
class StateSpace
{
public:

    virtual ~StateSpace();

    /**
     * Returns the size of the state space. 
     */
    virtual index_t size() const = 0;

    /**
     * Returns the coordinate system of the state space.
     */
    virtual const CoordinateSystem& crs() const = 0;
};



/**
 * Base class for state spaces with elements organized on a regular (multi-dimensional) grid.
 */
class GriddedStateSpace : public StateSpace
{
public:

    /**
     * Returns the number of grid dimensions. 
     * This is a convenient alias for `crs().ndim()`.
     */
    int ndim() const;


    /**
     * Returns the spatial extent of the grid.
     */
    virtual const AABox& extent() const = 0;

    /**
     * Returns the shape of the grid (the number of cells in each dimension).
     */
    virtual const ArrayShape& shape() const = 0;

    /**
     * Returns flat array of indexes identifying grid cells included in the state vector. 
     * 
     * If empty array is returned (`size()==0`), all grid cells are assumed to be included
     * in the state vector.
     */
    //virtual const Array& mask() const = 0;

};


/**
 * Abstract base class for state space partitioning schemes.
 * 
 * Implementations of this abstract base define how to partition the state space into local domains 
 * for localized analysis.
 */ 
class StateSpacePartitioning
{
public:

    StateSpacePartitioning(std::shared_ptr<const StateSpace> ss); 

    virtual ~StateSpacePartitioning();


    /** 
     * Returns the state space that is being partitioned.
     */
    virtual const StateSpace& stateSpace() const;

    /**
     * Returns the number of local analysis domains.
     */
    virtual int numDomains() const = 0;


    /** 
     * Returns the state vector size corresponding to the local domain `d`.
     */
    virtual index_t getLocalStateSize(int d) const = 0;


    /** 
     * Reads state vector for domain `d` from the global state.
     * 
     * @param d      Domain index.
     * @param Xg     Global state vector or ensemble.
     * @param outXl  Array where the local state vector or ensemble is to be stored.
     */
    virtual void getLocalState(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const = 0;


    /** 
     * Writes state vector for domain `d` into the global state.
     * 
     * @param d      Domain index.
     * @param Xl     Local state vector or ensemble for domain `d`.
     * @param Xg     Global state vector or ensemble.
     */
    virtual void putLocalState(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const = 0;


protected:
    std::shared_ptr<const StateSpace> mSS;
};



}

#endif