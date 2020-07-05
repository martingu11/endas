/**
 * @file GriddedPartitionPointQuery.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_GRIDDED_PARTITION_POINT_QUERY_HPP__
#define __ENDAS_DA_GRIDDED_PARTITION_POINT_QUERY_HPP__

#include <Endas/DA/DomainPartitioning.hpp>

#include <memory>

namespace endas
{

/** 
 * @addtogroup domain
 * @{ 
 */


/**
 * General purpose PartitionPointQuery implementation that indexes given points for faster 
 * retrieval.
 * 
 * The point query can work with any partitioning scheme that divides the domain into axis-aligned 
 * blocks and implements AABoxDomainPartitioning. The implementation internally uses kd-tree to 
 * organize the point data. 
 */ 
class IndexedPartitionPointQuery : public PartitionPointQuery
{
public:

    /**
     * IndexedPartitionPointQuery constructor.
     * 
     * Coordinates of the points to be indexed are given via the `coords` array, with point
     * coordinates stored in columns. Therefore, the first row correspomnds to the first coordinate, 
     * second row to the second coordinate etc. The number of columns defines the number of points
     * to index.
     * 
     * @param partitioner   Instance of the domain partitioner 
     * @param coords        Coordinates of the indexed point data
     */
    IndexedPartitionPointQuery(std::shared_ptr<const AABoxDomainPartitioning> partitioner,
                               Array2d coords); 
                            
    virtual ~IndexedPartitionPointQuery();

    virtual void rangeQuery(int domain, double range, IndexArray& out) const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};



/** @} */

}

#endif