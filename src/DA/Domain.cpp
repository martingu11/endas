#include <Endas/DA/Domain.hpp>
#include <Endas/DA/DomainPartitioning.hpp>

using namespace std;
using namespace endas;


DiscreteDomain::~DiscreteDomain()
{ }


int DiscreteSpatialDomain::dim() const
{
    return this->crs().dim();
}


bool GriddedDomain::hasEfficientSubset() const { return false; }


void GriddedDomain::getSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const
{
    IndexArray indices;
    this->getIndices(block, indices);
    selectRows(X, indices, out);
}

void GriddedDomain::putSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const
{
    IndexArray indices;
    this->getIndices(block, indices);
    distributeRows(X, indices, out);
}


DomainPartitioning::~DomainPartitioning()
{ }