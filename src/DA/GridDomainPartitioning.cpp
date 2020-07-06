#include <Endas/DA/GridDomainPartitioning.hpp>
#include <Endas/DA/IndexedPartitionPointQuery.hpp>
#include <Endas/Endas.hpp>

#include "../Compatibility.hpp"

#include <Eigen/Geometry>

using namespace std;
using namespace endas;


struct Domain
{
    GriddedDomain::Block block;
    IndexArray indices;
    index_t size;

    Domain(index_t dim): block(dim) { }
};







struct GridDomainPartitioning::Data
{
    shared_ptr<const GriddedDomain> globalDomain;
    int blockSize;
    int padding;

    // Local domains.
    vector<Domain> domains;

    void generateDomains();
};


GridDomainPartitioning::GridDomainPartitioning(shared_ptr<const GriddedDomain> domain, 
                                               int blockSize, int padding) 
: mData(make_unique<Data>())
{ 
    ENDAS_ASSERT(domain);
    ENDAS_REQUIRE(domain->coordDim() > 0 && domain->coordDim() <= 2, 
        std::invalid_argument, "Only one or two-dimensional grids can currently be partitioned");

    mData->globalDomain = domain;
    mData->blockSize = blockSize;
    mData->padding = padding;

    mData->generateDomains();
}

GridDomainPartitioning::~GridDomainPartitioning()
{ }


const DiscreteDomain& GridDomainPartitioning::domain() const
{
    ENDAS_ASSERT(mData->globalDomain);
    return *mData->globalDomain;
}


int GridDomainPartitioning::partitionCoordDim() const
{
    ENDAS_ASSERT(mData->globalDomain);
    return mData->globalDomain->coordDim();
}


void GridDomainPartitioning::Data::generateDomains()
{
    ENDAS_ASSERT(globalDomain);
    int dim = globalDomain->coordDim();
    int useIndices = !globalDomain->hasEfficientSubset();

    auto addDomainFn = [&](Domain&& domain)
    {
        if (useIndices) 
        {
            globalDomain->getIndices(domain.block, domain.indices);
            domain.size = domain.indices.size();
        }
        else
            domain.size = globalDomain->blockSize(domain.block);

        if (domain.size > 0) domains.push_back(domain);
    };


    // One-dimensional case
    if (dim == 1)
    {
        index_t N = globalDomain->shape()(0);
        for (index_t i = 0; i < N; i+= blockSize)
        {
            Domain domain(1);
            domain.block.min()(0) = i;
            domain.block.max()(0) = std::min(i+blockSize, N);
            addDomainFn(move(domain));
        }
    }
    else if (dim == 2)
    {
        index_t NX = globalDomain->shape()(0);
        index_t NY = globalDomain->shape()(1);

        for (index_t x = 0; x < NX; x+= blockSize)
        {
            index_t xend = std::min(x+blockSize, NX);
            for (index_t y = 0; y < NY; y+= blockSize)
            {
                Domain domain(2);
                domain.block.min()(0) = x;
                domain.block.min()(1) = y;
                domain.block.max()(0) = xend;
                domain.block.max()(1) = std::min(y+blockSize, NY);
                addDomainFn(move(domain));
            }
        }
    }
    else
    {
        ENDAS_ASSERT(false);
    }
}


int GridDomainPartitioning::numLocalDomains() const
{
    return (int)mData->domains.size();
}

index_t GridDomainPartitioning::getLocalSize(int d) const
{
    ENDAS_ASSERT(d >= 0 && d < mData->domains.size());
    return mData->domains[d].size;
}

AABox GridDomainPartitioning::getLocalBox(int d) const
{
    ENDAS_ASSERT(d >= 0 && d < mData->domains.size());
    return mData->globalDomain->getBlockExtent(mData->domains[d].block);
}


void GridDomainPartitioning::getLocal(int d, const Ref<const Array2d> Xg, Ref<Array2d> out) const
{
    ENDAS_ASSERT(d >= 0 && d < mData->domains.size());
    const auto& domain = mData->domains[d];

    // Empty domains are allowed
    if (domain.size == 0) return;

    // Using indices
    if (domain.indices.size() > 0)
    {
        selectRows(Xg, domain.indices, out);
    }
    // Using getSubset()
    else
    {
        mData->globalDomain->getSubset(domain.block, Xg, out);
    }
}

void GridDomainPartitioning::putLocal(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const
{
    ENDAS_ASSERT(d >= 0 && d < mData->domains.size());
    const auto& domain = mData->domains[d];

    // Empty domains are allowed
    if (domain.size == 0) return;


    // Using indices
    if (domain.indices.size() > 0)
    {
        distributeRows(Xl, domain.indices, Xg);
    }
    // Using getSubset()
    else
    {
        mData->globalDomain->putSubset(domain.block, Xl, Xg);
    }

}


shared_ptr<const PartitionPointQuery> GridDomainPartitioning::indexPoints(Array2d coords) const
{
    return make_shared<IndexedPartitionPointQuery>(shared_ptr_wrap(*this), move(coords));
}

