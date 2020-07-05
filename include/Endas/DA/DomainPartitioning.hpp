/**
 * @file DomainPartitioning.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_DOMAIN_PARTITIONING_HPP__
#define __ENDAS_DA_DOMAIN_PARTITIONING_HPP__

#include "Domain.hpp"
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
 * domains of a partitioned domain. PartitionPointQuery instances are usually not created 
 * directly but are obtained from the domain partitioner. Their intended use is in observation
 * manager implementations where they can help with locating observations needed for a local domain.
 *  
 * @rst
 * .. tip::
 *    If the default querying capabilities provided by the domain partitioner do not fit your 
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
 * Abstract base class for domain partitioning schemes.
 * 
 * Implementations of this abstract base define how to partition a domain into local domains 
 * for localized analysis.
 */ 
class ENDAS_DLL DomainPartitioning
{
public:

    virtual ~DomainPartitioning();

    /** 
     * Returns the space that is being partitioned.
     */
    virtual const DiscreteDomain& domain() const = 0;


    /**
     * Returns the number of *partition* coordinate dimensions.
     * 
     * The partition dimension refers to how many coordinates are necessary to describe any location 
     * in the subspace of the domain in which the partitioning occurs. As an example, for a state 
     * space organized on two-dimensional grid and partitioned along the two spatial dimensions, we 
     * would use two coordinates to represent both the location of state space elements as well as
     * boundaries of the partitions. One could however have three-dimensional domain partitioned in
     * two dimensions only, in this case the coordinate dimension of the domain and of the partitioner
     * are three and two, respenctively.
     * 
     * See the documentation of individual domain partitioning schemes for more information.
     */
    virtual int partitionCoordDim() const = 0;


    /**
     * Returns the number of local analysis domains.
     */
    virtual int numLocalDomains() const = 0;


    /** 
     * Returns the state vector size corresponding to the local domain `d`.
     */
    virtual index_t getLocalSize(int d) const = 0;


    /** 
     * Reads state vector for local domain `d` from the global domain.
     * 
     * @param d      Domain index.
     * @param Xg     Global state vector or ensemble.
     * @param outXl  Array where the local state vector or ensemble is to be stored.
     */
    virtual void getLocal(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const = 0;


    /** 
     * Writes state vector for local domain `d` into the global domain.
     * 
     * @param d      Domain index.
     * @param Xl     Local state vector or ensemble for domain `d`.
     * @param Xg     Global state vector or ensemble.
     */
    virtual void putLocal(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const = 0;


    /**
     * Indexes given point data and returns query implementation.
     * 
     * The returned PartitionPointQuery instace can be used to retrieve subsets of the indexed
     * points for individual local domains generated by this partitioning scheme. The indexed 
     * points are given as a two-dimensional array and are stored in columns, i.e. the array 
     * should have one row per coordinate dimension (should match partitionCoordDim()) and one column per
     * indexed point.
     * 
     * @param coords Point coordinates to be indexed. 
     */
    virtual std::shared_ptr<const PartitionPointQuery> indexPoints(Array2d coords) const = 0;

};


/**
 * Abstract base class for domain partitioning schemes based on bounding boxes.
 * 
 * AABoxDomainPartitioning represents partitioning schemes that divide the global domain into local 
 * domains that are spatially bounded by an axis-aligned box.
 */
class ENDAS_DLL AABoxDomainPartitioning : public DomainPartitioning
{
public:

    /** 
     * Returns bounding box of the local domain `d`.
     */
    virtual AABox getLocalBox(int d) const = 0;

};


/** @} */

}

#endif