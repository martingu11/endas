#include <Endas/DA/Domain.hpp>
#include <Endas/DA/DomainPartitioning.hpp>

using namespace std;
using namespace endas;


DiscreteDomain::~DiscreteDomain()
{ }


int DiscreteSpatialDomain::coordDim() const
{
    return this->crs().dim();
}

void DiscreteSpatialDomain::getCoords(const IndexArray& selected, Ref<Array2d> out) const
{
    auto dim = this->coordDim();
    ENDAS_ASSERT(out.rows() == dim);
    ENDAS_ASSERT(out.cols() == selected.size());

    Array2d allCoords(dim, this->size());
    this->getCoords(allCoords);

    selectCols(allCoords, selected, out);
}


Array GriddedDomain::cellSize() const
{
    return (this->extent().max() - this->extent().min()).array() / this->shape().cast<AABox::Scalar>();
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