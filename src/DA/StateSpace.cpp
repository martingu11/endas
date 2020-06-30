#include <Endas/DA/StateSpacePartitioning.hpp>

using namespace std;
using namespace endas;


StateSpace::~StateSpace()
{ }

int SpatialStateSpace::dim() const
{
    return this->crs().dim();
}


bool GriddedStateSpace::hasEfficientSubset() const { return false; }


void GriddedStateSpace::getSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const
{
    IndexArray indices;
    this->getIndices(block, indices);
    selectRows(X, indices, out);
}

void GriddedStateSpace::putSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const
{
    IndexArray indices;
    this->getIndices(block, indices);
    distributeRows(X, indices, out);
}




StateSpacePartitioning::~StateSpacePartitioning()
{ }