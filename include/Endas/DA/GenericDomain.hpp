/**
 * @file GenericDomain.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DOMAIN_GENERIC_DOMAIN_HPP__
#define __ENDAS_DOMAIN_GENERIC_DOMAIN_HPP__

#include <Endas/DA/Domain.hpp>
#include <Endas/DA/DomainPartitioning.hpp>
#include <Endas/Spatial/CoordinateSystem.hpp>

namespace endas
{

/** 
 * @addtogroup domain
 * @{ 
 */


/**
 * Trivial domain without assumptions on the structure.
 *
 * GenericDomain is likely only useful for low-dimensional toy problems where each state variable
 * can be updated independently. Both the DiscreteDomain and DomainPartitioning interfaces are 
 * implemented by this class as there is very little benefit in separating the two.
 * 
 * Under the implemented partitioning scheme, the index of each domain element is also its coordinate
 * and coordDim() is therefore 1. One-dimensional Euclidean distance is used by default to evaluate 
 * proximity although this can be replaced by user-defined function.
 */ 
class GenericDomain : public DiscreteDomain, public DomainPartitioning
{
public:

    GenericDomain(index_t size);

    virtual index_t size() const override;
    virtual int coordDim() const override;    
    virtual const DiscreteDomain& domain() const override;

    virtual int numLocalDomains() const override;
    virtual index_t getLocalSize(int d) const override;
    virtual void getLocal(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const override;
    virtual void putLocal(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const override;

    virtual std::shared_ptr<const PartitionPointQuery> 
    indexPoints(const Ref<const Array2d> coords) const override;

private:
    index_t mSize;
};



/** @} */

}

#endif