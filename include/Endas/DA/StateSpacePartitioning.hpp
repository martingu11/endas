/**
 * @file StateSpacePartitioning.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_STATESPACE_PARTITIONING_HPP__
#define __ENDAS_DA_STATESPACE_PARTITIONING_HPP__

#include "StateSpace.hpp"
#include <vector>

namespace endas
{

/** 
 * @addtogroup da
 * @{ 
 */


/**
 * Efficient querying of observations for localized analysis.
 * 
 * PartitionPointQuery is an abstract query mechanism for retrieving observations for local analysis
 * domains of a partitioned state space. PartitionPointQuery instances are usually not created 
 * directly but are obtained from the state space partitioner. Their intended use is in observation
 * manager implementations where they can help with locating observations needed for a local domain.
 *  
 * @rst
 * .. tip::
 *    If the default querying capabilities provided by the state space partitioner do not fit your 
 *    needs, you can roll out your completely own and new observation manager. This gives you the 
 *    freedom to follow any approach to locate and fetch observations.
 * @endrst
 */
class ENDAS_DLL PartitionPointQuery
{
public:
  
    virtual ~PartitionPointQuery() { }

    /**
     * Queries points within the given distance from a local domain.
     * 
     * @param d      Local domain index
     * @param range  The search range
     * @param out    Array of indices to populate
     */ 
    virtual void rangeQuery(int domain, double range, IndexArray& out) const = 0;

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

    virtual ~StateSpacePartitioning();

    /** 
     * Returns the state space that is being partitioned.
     */
    virtual const StateSpace& stateSpace() const = 0;


    /**
     * Returns the number of *partition* coordinate dimensions.
     * 
     * The partition dimension refers to how many coordinates are necessary to describe any location 
     * in the subspace of the state space in which the partitioning occurs. As an example, for a state 
     * space organized on two-dimensional grid and partitioned along the two spatial dimensions, we 
     * would most likely use two coordinates to represent any location within the partitioning subspace.
     * One could however choose more esoteric coordinate representation as well.
     * 
     * See the documentation of individual state space partitioning schemes for more information.
     */
    virtual int coordDim() const = 0;


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


    /**
     * Indexes given point data and returns query implementation.
     * 
     * The returned PartitionPointQuery instace can be used to retrieve subsets of the indexed
     * points for individual local domains generated by this partitioning scheme. The indexed 
     * points are given as a two-dimensional array and are stored in columns, i.e. the array 
     * should have one row per coordinate dimension (should match coordDim()) and one column per
     * indexed point.
     * 
     * @param coords Point coordinates to be indexed. 
     */
    virtual std::shared_ptr<const PartitionPointQuery> indexPoints(const Ref<const Array2d> coords) const = 0;

};


/** @} */

}

#endif