#include <Endas/DA/GridDomain.hpp>
#include "../Compatibility.hpp"

#include <iostream>

using namespace std;
using namespace endas;


GridDomain::GridDomain(const ArrayShape& shape, std::shared_ptr<const CoordinateSystem> crs, 
                       const AABox& extent, int numVarsPerCell)
: mShape(shape), mCRS(crs), mExtent(extent), mNumVarsPerCell(numVarsPerCell)
{ 
    ENDAS_ASSERT(shape.size() == crs->dim());
    ENDAS_ASSERT(numVarsPerCell > 0);

    mSize = shape.prod() * numVarsPerCell;
    mCellSize = (extent.max() - extent.min()).array() / shape.cast<AABox::Scalar>();
}


GridDomain::GridDomain(const ArrayShape& shape, std::shared_ptr<const CoordinateSystem> crs, 
                       const AABox& extent, IndexArray cellMap)
: mShape(shape), mCRS(crs), mExtent(extent), mNumVarsPerCell(0), mCellMap(move(cellMap))
{
    ENDAS_ASSERT(shape.size() == crs->dim());
    ENDAS_ASSERT(mCellMap.size() > 0);

    mSize = mCellMap.size();
    mCellSize = (extent.max() - extent.min()).array() / shape.cast<AABox::Scalar>();
}


GridDomain::~GridDomain()
{ }

index_t GridDomain::size() const
{
    return mSize;
}

const ArrayShape& GridDomain::shape() const
{
    return mShape;
}


const CoordinateSystem& GridDomain::crs() const
{
    ENDAS_ASSERT(mCRS);
    return *mCRS;
}

const AABox& GridDomain::extent() const
{
    return mExtent;
}


const IndexArray& GridDomain::cellMap() const
{
    return mCellMap;
}


index_t GridDomain::blockSize(const GriddedDomain::Block& block) const
{
    ENDAS_ASSERT(block.min().minCoeff() >= 0);
    ENDAS_ASSERT((block.max().array() <= mShape).all());

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

AABox GridDomain::getBlockExtent(const Block& block) const
{
    auto min = mExtent.min().array() + (mCellSize * block.min().array().cast<real_t>());
    auto max = mExtent.min().array() + (mCellSize * block.max().array().cast<real_t>());
    return AABox(min, max);
}



// Calls fn(i, iend) for each continuous range of state variables within a block, where `i` and 
// `iend` denote the continuous range of state variables that are indide the block.
template <class Fn> inline void 
forEachBlockStateRange(const GriddedDomain::Block& block, int numVarsPerCell, 
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
        index_t ysize = gridShape(1) * numVarsPerCell;

        for (int x = block.min()(0); x != block.max()(0); x++) 
        {
            index_t i = x * ysize + block.min()(1) * numVarsPerCell;
            index_t iend = x * ysize + block.max()(1) * numVarsPerCell;
            fn(i, iend);
        }
    }
    else 
    {
        ENDAS_NOT_IMPLEMENTED;
    }
}

void GridDomain::getIndices(const Block& block, IndexArray& out) const
{
    // Dense grid
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



void GridDomain::getCoords(Ref<Array2d> out) const
{
    int n = this->size();
    int dim = this->coordDim();
    ENDAS_ASSERT(out.rows() >= dim); 
    ENDAS_ASSERT(out.cols() >= n); 

    // Dense grid
    if (mCellMap.size() == 0)
    {
        Block block(Block::VectorType::Zero(dim), mShape);
        forEachBlockStateRange(block, mNumVarsPerCell, mShape, [&](index_t i, index_t iend)
        {
            while (i != iend) 
            {
                cellCoord(i / mNumVarsPerCell, out.col(i));
                i++;
            }
        });
    }
    else
    {
        ENDAS_NOT_IMPLEMENTED;
    }
}

void GridDomain::getCoords(const IndexArray& selected, Ref<Array2d> out) const
{
    int n = selected.size();
    int dim = this->coordDim();
    ENDAS_ASSERT(out.rows() >= dim); 
    ENDAS_ASSERT(out.cols() >= n); 

    // Dense grid
    if (mCellMap.size() == 0)
    {
        for (int i = 0; i != n; i++)
        {
            cellCoord(selected[i] / mNumVarsPerCell, out.col(i));
        }
    }
    else
    {
        ENDAS_NOT_IMPLEMENTED;
    }
}


void GridDomain::cellCoord(index_t cellIndex, Ref<Array> out) const
{
    auto dim = this->coordDim();
    
    if (dim == 1)
    {
        out(0) = mExtent.min()(0) + mCellSize(0) * cellIndex;
    }
    else if (dim == 2)
    {
        index_t xi = (index_t)(cellIndex / mShape(0));
        index_t yi = (index_t)(cellIndex % mShape(0));
        out(0) = mExtent.min()(0) + mCellSize(0)*xi;        
        out(1) = mExtent.min()(1) + mCellSize(1)*yi;
    }
    else 
    {
        ENDAS_NOT_IMPLEMENTED;
    }

}




bool GridDomain::hasEfficientSubset() const
{
    return mCellMap.size() == 0; // Not using cellmap
}


void GridDomain::getSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const
{
    ENDAS_ASSERT(X.cols() == out.cols());

    if (mCellMap.size() == 0)
    {
        index_t ilocal = 0;
        forEachBlockStateRange(block, mNumVarsPerCell, mShape, [&](index_t i, index_t iend)
        {
            index_t n = iend - i;
            index_t N = X.cols();

            /*cout << "    getsubset " 
                << ilocal << " 0 " << n << " " << N << " <- " 
                << i << " 0 " << n << " " << N << endl;*/

            out.block(ilocal, 0, n, N) = X.block(i, 0, n, N);
            ilocal += n;
        });
    }
    else
    {
        GriddedDomain::getSubset(block, X, out);
    }
}


void GridDomain::putSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const
{
    ENDAS_ASSERT(X.cols() == out.cols());

    if (mCellMap.size() == 0)
    {
        index_t ilocal = 0;
        forEachBlockStateRange(block, mNumVarsPerCell, mShape, [&](index_t i, index_t iend)
        {
            index_t n = iend - i;
            index_t N = X.cols();

            /*cout << "    putsubset " 
                << ilocal << " 0 " << n << " " << N << " -> " 
                << i << " 0 " << n << " " << N << endl;*/

            out.block(i, 0, n, N) = X.block(ilocal, 0, n, N);
            ilocal+= n;
        });
    }
    else
    {
        GriddedDomain::putSubset(block, X, out);
    }
}



