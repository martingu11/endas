#include <Endas/Domain/GridStateSpacePartitioning.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"


#include <Eigen/Geometry>

using namespace std;
using namespace endas;


struct Domain
{
    GriddedStateSpace::Block block;
    IndexArray indices;
    index_t size;

    Domain(index_t dim): block(dim) { }
};

struct GridStateSpacePartitioning::Data
{
    shared_ptr<const GriddedStateSpace> ss;
    int blockSize;
    int padding;

    // State space partitions/domains.
    vector<Domain> domains;

    void generateDomains();


};


GridStateSpacePartitioning::GridStateSpacePartitioning(shared_ptr<const GriddedStateSpace> stateSpace, 
                                                       int blockSize, int padding) 
: mData(make_unique<Data>())
{ 
    ENDAS_REQUIRE(stateSpace->dim() > 0 && stateSpace->dim() <= 2, 
        std::invalid_argument, "Only one or two-dimensional grids can be partitioned");

    mData->ss = stateSpace;
    mData->blockSize = blockSize;
    mData->padding = padding;

    mData->generateDomains();
}

GridStateSpacePartitioning::~GridStateSpacePartitioning()
{ }


void GridStateSpacePartitioning::Data::generateDomains()
{
    int dim = ss->dim();
    int useIndices = !ss->hasEfficientSubset();

    auto addDomainFn = [&](Domain&& domain)
    {
        if (useIndices) 
        {
            ss->getIndices(domain.block, domain.indices);
            domain.size = domain.indices.size();
        }
        else
            domain.size = ss->size(domain.block);

        if (domain.size > 0) domains.push_back(domain);
    };


    // One-dimensional case
    if (dim == 1)
    {
        index_t N = ss->shape()(0);
        for (index_t i = 0; i != N; i+= blockSize)
        {
            Domain domain(1);
            domain.block.min()(0) = i;
            domain.block.max()(0) = std::min(i+blockSize, N);
            addDomainFn(move(domain));
        }
    }
    else if (dim == 2)
    {
        index_t NR = ss->shape()(0);
        index_t NC = ss->shape()(1);

        for (index_t i = 0; i != NR; i+= blockSize)
        {
            index_t iend = std::min(i+blockSize, NR);
            for (index_t j = 0; j != NC; j+= blockSize)
            {
                Domain domain(2);
                domain.block.min()(0) = i;
                domain.block.min()(1) = j;
                domain.block.max()(0) = iend;
                domain.block.max()(1) = std::min(j+blockSize, NC);
                addDomainFn(move(domain));
            }
        }
    }
    else
    {
        ENDAS_ASSERT(false);
    }
}


int GridStateSpacePartitioning::numDomains() const
{
    return (int)mData->domains.size();
}

index_t GridStateSpacePartitioning::getLocalStateSize(int d) const
{
    ENDAS_ASSERT(d >= 0 && d < mData->domains.size());
    return mData->domains[d].size;
}

void GridStateSpacePartitioning::getLocalState(int d, const Ref<const Array2d> Xg, 
                                               Ref<Array2d> out) const
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
        mData->ss->getSubset(domain.block, Xg, out);
    }
}

void GridStateSpacePartitioning::putLocalState(int d, const Ref<const Array2d> Xl, 
                                               Ref<Array2d> Xg) const
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
        mData->ss->putSubset(domain.block, Xl, Xg);
    }

}


shared_ptr<const PartitionPointQuery> 
GridStateSpacePartitioning::indexPoints(const Ref<const Array2d> coords) const
{
    ENDAS_NOT_IMPLEMENTED;
}

