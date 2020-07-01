/**
 * @file GridDomainPartitioning.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_GRID_DOMAIN_PARTITIONING_HPP__
#define __ENDAS_DA_GRID_DOMAIN_PARTITIONING_HPP__

#include <Endas/DA/Domain.hpp>
#include <Endas/DA/DomainPartitioning.hpp>

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
class GridDomainPartitioning : public DomainPartitioning
{
public:

    GridDomainPartitioning(std::shared_ptr<const GriddedDomain> stateSpace, 
                           int blockSize, int padding = 0); 
                            
    virtual ~GridDomainPartitioning();

    virtual int numLocalDomains() const override;
    virtual index_t getLocalSize(int d) const override;
    virtual void getLocal(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const override;
    virtual void putLocal(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const override;

    virtual std::shared_ptr<const PartitionPointQuery> 
    indexPoints(const Ref<const Array2d> coords) const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};



/** @} */

}

#endif