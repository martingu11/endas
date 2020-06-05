/**
 * @file StateSpaceGrid.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_GRID_STATE_SPACE_PARTITIONING_HPP__
#define __ENDAS_DA_GRID_STATE_SPACE_PARTITIONING_HPP__

#include <Endas/DA/StateSpace.hpp>
#include <Endas/DA/StateSpacePartitioning.hpp>

#include <memory>

namespace endas
{

/** 
 * @addtogroup domain
 * @{ 
 */


/**
 * State space partitioning scheme that operates on gridded state spaces.
 * 
 * The partitioning scheme divides the state space into rectangular local domains of fixed size.
 * Domains may be fully disjoint or partly overlap. In the latter case, the overlapping regions of
 * adjacent local domains are blended together to remove any visible boundaries in the global
 * state after observations have been assimilated.
 * 
 */ 
class GridStateSpacePartitioning : public StateSpacePartitioning
{
public:

    GridStateSpacePartitioning(std::shared_ptr<const GriddedStateSpace> stateSpace, 
                               int blockSize, int padding = 0); 
                            
    virtual ~GridStateSpacePartitioning();

    virtual int numDomains() const override;
    virtual index_t getLocalStateSize(int d) const override;
    virtual void getLocalState(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const override;
    virtual void putLocalState(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const override;

    virtual std::shared_ptr<const PartitionPointQuery> 
    indexPoints(const Ref<const Array2d> coords) const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};



/** @} */

}

#endif