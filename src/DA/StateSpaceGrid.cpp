#include <Endas/DA/StateSpaceGrid.hpp>
#include "../Compatibility.hpp"


using namespace std;
using namespace endas;


struct StateSpaceGrid::Data
{
    ArrayShape shape;
    AABox extent;
    shared_ptr<const CoordinateSystem> crs;

    index_t size;
};


StateSpaceGrid::StateSpaceGrid(const ArrayShape& shape, std::shared_ptr<const CoordinateSystem> crs, 
                               const AABox& extent)
: mData(make_unique<Data>())
{ 
    if (shape.size() > 2)
    {
        ENDAS_NOT_SUPPORTED("StateSpaceGrid only supports up to 2 spatial dimensions");
    }

    ENDAS_ASSERT(shape.size() == crs->ndim());

    mData->shape = shape;
    mData->extent = extent;
    mData->crs = crs;
    mData->size = shape.prod();
}


StateSpaceGrid::~StateSpaceGrid()
{ }

index_t StateSpaceGrid::size() const
{
    return mData->size;
}

const ArrayShape& StateSpaceGrid::shape() const
{
    return mData->shape;
}


const CoordinateSystem& StateSpaceGrid::crs() const
{
    ENDAS_ASSERT(mData->crs);
    return *mData->crs;
}

const AABox& StateSpaceGrid::extent() const
{
    return mData->extent;
}


