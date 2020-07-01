/**
 * @file GridDomain.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_GRID_DOMAIN_HPP__
#define __ENDAS_DA_GRID_DOMAIN_HPP__

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
 * Generic implementation of a discrete domain with elements organized on a multi-dimensional grid.
 * 
 * The grid can be either dense with fixed number of variables per cell, or state variables can
 * be assigned to grid cells arbitrarily via a provided mapping. For dense grids, the state vector 
 * must be organized so that grid cells are laid out in column-major order in the vector (one dimensional 
 * dense grids are assumed to have a single column). Furthermore, variables for a single cell must 
 * be adjacent in the state vector and always in the same order. When the state variable-to-grid-cell
 * mapping is provided, no assumptions are made on how the state vector is organized. Please note 
 * that ``hasEfficientSubset()`` will return ``false`` for mapped grids.
 * 
 *  
 * @rst
 * .. figure:: /images/loc_grid2d_masking.*
 * 
 *    Dense and mapped variants demonstrated on a two-dimensional grid.
 *
 * .. note:: 
 *    One and two -dimensional grids are currently supported when the dense option is in use. 
 *    This limitation does not exist for mapped grids.
 * @endrst
 */ 
class GridDomain : public GriddedDomain
{
public:


     /**
     * Creates dense StateSpaceGrid.
     * 
     * The grid will have the same number of state variables in each cells and the state vector 
     * must be organized so that grid cells are laid out column-by-column. Furthermore, variables 
     * for a single cell must be adjacent in the state vector and always in the same order.
     * 
     * @param shape          The size of the grid in each dimension.
     * @param crs            Coordinate system of the grid. The instance is copied.
     * @param extent         Physical extent (bounding box) of the grid.
     * @param numVarsPerCell The number of state variables in each grid cell.
     */
    GridDomain(const ArrayShape& shape, std::shared_ptr<const CoordinateSystem> crs, 
                   const AABox& extent, int numVarsPerCell = 1);

   /**
     * Creates dense StateSpaceGrid.
     * 
     * This is a convenience overload that copies the CRS. 
     * 
     * See ``GridDomain::GridDomain(shape, crs, extent, numVarsPerCell)`` for explanation
     * of the parameters.
     */
    template <class CRS, require_is_convertible<CRS, const CoordinateSystem> = true>
    GridDomain(const ArrayShape& shape, const CRS& crs, const AABox& extent, int numVarsPerCell = 1)
    : GridDomain(shape, std::make_shared<CRS>(crs), extent, numVarsPerCell)
    { }


    /**
     * Creates StateSpaceGrid with arbitrary mapping to the state vector.
     * 
     * The grid is constructed from a mapping that assigns each state variable to a single grid cell.
     * No assumptions are made on how the state vector is organized. The ``cellMap`` array must have 
     * the same size as the state vector and each element must hold the index of the corresponding 
     * grid cell. Cells are indexed from top-to-bottom and left-to-right (i.e. column-wise).
     * 
     * @rst
     * .. tip::
     *    You can use ``std::move`` to avoid copying the cell map::
     * 
     *        IndexArray cellMap = { ... };  
     *        GridDomain gss(shape, crs, extent, std::move(cellMap));
     *        // At this point cellMap is invalid, use gss.cellMap() instead
     * @endrst
     * 
     * 
     * @param shape      The size of the grid in each dimension.
     * @param crs        Coordinate system of the grid. 
     * @param extent     Physical extent (bounding box) of the grid.
     * @param cellMap    Index array mapping each state variable to a grid cell.
     */
    GridDomain(const ArrayShape& shape, std::shared_ptr<const CoordinateSystem> crs, 
                   const AABox& extent, IndexArray cellMap);


    /**
     * Creates StateSpaceGrid with arbitrary mapping to the state vector.
     * 
     * This is a convenience overload that copies the CRS. 
     * 
     * See ``GridDomain::GridDomain(shape, crs, extent, cellMap)`` for explanation of the 
     * parameters.
     */
    template <class CRS, require_is_convertible<CRS, const CoordinateSystem> = true>
    GridDomain(const ArrayShape& shape, const CRS& crs, const AABox& extent, IndexArray cellMap)
    : GridDomain(shape, std::make_shared<CRS>(crs), extent, cellMap)
    { }


    virtual ~GridDomain();

    virtual index_t size() const override;
    virtual const CoordinateSystem& crs() const override;
    virtual const AABox& extent() const override;
    virtual const ArrayShape& shape() const override;

    /** Returns the cell map or empty array if not used. */
    const IndexArray& cellMap() const;

    virtual index_t size(const GriddedDomain::Block& block) const override;

    virtual void getCoords(Ref<Array2d> out) const override;

    virtual void getIndices(const Block& block, IndexArray& out) const override;
    virtual bool hasEfficientSubset() const override;
    virtual void getSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const override;
    virtual void putSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const override;



private:
    index_t mSize;
    ArrayShape mShape;
    AABox mExtent;
    Array mCellSize;

    std::shared_ptr<const CoordinateSystem> mCRS;

    int mNumVarsPerCell;      // For dense grids
    IndexArray mCellMap;      // For arbitrary grids


    void cellCoord(index_t cellIndex, Ref<Array> out) const;

};


/** @} */

}

#endif