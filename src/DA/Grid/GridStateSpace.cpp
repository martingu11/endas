#include <Endas/DA/Grid/GridStateSpace.hpp>
#include "../../Compatibility.hpp"


using namespace std;
using namespace endas;


GridStateSpace::GridStateSpace(const ArrayShape& shape, std::shared_ptr<const CoordinateSystem> crs, 
                               const AABox& extent, int numVarsPerCell)
: mShape(shape), mCRS(crs), mExtent(extent), mNumVarsPerCell(numVarsPerCell)
{ 
    ENDAS_ASSERT(shape.size() == crs->dim());
    ENDAS_ASSERT(numVarsPerCell > 0);

    mSize = shape.prod() * numVarsPerCell;
}


GridStateSpace::GridStateSpace(const ArrayShape& shape, std::shared_ptr<const CoordinateSystem> crs, 
                               const AABox& extent, IndexArray cellMap)
: mShape(shape), mCRS(crs), mExtent(extent), mNumVarsPerCell(0), mCellMap(move(cellMap))
{
    ENDAS_ASSERT(shape.size() == crs->dim());
    ENDAS_ASSERT(mCellMap.size() > 0);

    mSize = mCellMap.size();
}


GridStateSpace::~GridStateSpace()
{ }

index_t GridStateSpace::size() const
{
    return mSize;
}

const ArrayShape& GridStateSpace::shape() const
{
    return mShape;
}


const CoordinateSystem& GridStateSpace::crs() const
{
    ENDAS_ASSERT(mCRS);
    return *mCRS;
}

const AABox& GridStateSpace::extent() const
{
    return mExtent;
}


const IndexArray& GridStateSpace::cellMap() const
{
    return mCellMap;
}


index_t GridStateSpace::size(const GriddedStateSpace::Block& block) const
{
    ENDAS_ASSERT(block.min().minCoeff() > 0);
    ENDAS_ASSERT((block.max().array().cast<int>() <= mShape).all());

    // Dense grid
    if (mCellMap.size() == 0)
    {
        return block.volume() * mNumVarsPerCell;
    }
    else
    {
        ENDAS_NOT_IMPLEMENTED;
    }
}


// Calls fn(i, iend) for each continuous range of state variables within a block, where `i` and 
// `iend` denote the continuous range of state variables that are indide the block.
template <class Fn> inline void 
forEachBlockStateRange(const GriddedStateSpace::Block& block, int numVarsPerCell, 
                       const ArrayShape& gridShape, Fn fn)
{
    int dim = gridShape.size();
    ENDAS_ASSERT(dim == block.dim());

    if (dim == 1)
    {
        index_t i = block.min()(0) * numVarsPerCell;
        index_t iend = block.max()(0) * numVarsPerCell;
        fn(i, iend);
    }
    else if (dim == 2)
    {
        index_t colsize = gridShape(0) * numVarsPerCell;
        for (int j = block.min()(1); j != block.max()(1); j++) 
        {
            index_t i = j * colsize + block.min()(0) * numVarsPerCell;
            index_t iend = j * colsize + block.max()(0) * numVarsPerCell;
            fn(i, iend);
        }
    }
    else 
    {
        ENDAS_NOT_IMPLEMENTED;
    }
}

void GridStateSpace::getIndices(const Block& block, IndexArray& out) const
{
    // Not using cellmap
    if (mCellMap.size() == 0) 
    {
        forEachBlockStateRange(block, mNumVarsPerCell, mShape, [&](index_t i, index_t iend)
        {
            while (i != iend) out.push_back(i++);
        });
    }
    else 
    {
        ENDAS_NOT_IMPLEMENTED;
    }
}


bool GridStateSpace::hasEfficientSubset() const
{
    return mCellMap.size() == 0; // Not using cellmap
}


void GridStateSpace::getSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const
{
    ENDAS_ASSERT(X.cols() == out.cols());

    if (mCellMap.size() == 0)
    {
        forEachBlockStateRange(block, mNumVarsPerCell, mShape, [&](index_t i, index_t iend)
        {
            out = X.block(i, 0, iend, X.cols());
        });
    }
    else
    {
        GriddedStateSpace::getSubset(block, X, out);
    }
}


void GridStateSpace::putSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const
{
    ENDAS_ASSERT(X.cols() == out.cols());

    if (mCellMap.size() == 0)
    {
        forEachBlockStateRange(block, mNumVarsPerCell, mShape, [&](index_t i, index_t iend)
        {
            out.block(i, 0, iend, X.cols()) = X;
        });
    }
    else
    {
        GriddedStateSpace::putSubset(block, X, out);
    }
}



