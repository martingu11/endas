/**
 * @file StateSpaceGrid.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_STATESPACE_GRID_HPP__
#define __ENDAS_DA_STATESPACE_GRID_HPP__

#include "StateSpace.hpp"
#include <memory>

namespace endas
{


/**
 * State space with elements organized on a multi-dimensional grid.
 * 
 */ 
class StateSpaceGrid : public GriddedStateSpace
{
public:


    /**
     * StateSpaceGrid constructor.
     * 
     * @param shape   The size of the grid in each dimension.
     * @param crs     Coordinate system of the grid. The instance is copied.
     * @param extent  Physical extent (bounding box) of the grid.
     */
    template <class CRS, require_is_convertible<CRS, const CoordinateSystem> = true>
    StateSpaceGrid(const ArrayShape& shape, const CRS& crs, const AABox& extent)
    : StateSpaceGrid(shape, std::make_shared<CRS>(crs), extent)
    { }


    /**
     * StateSpaceGrid constructor.
     * 
     * @param shape   The size of the grid in each dimension.
     * @param crs     Coordinate system of the grid.
     * @param extent  Physical extent (bounding box) of the grid.
     */
    StateSpaceGrid(const ArrayShape& shape, std::shared_ptr<const CoordinateSystem> crs, 
                   const AABox& extent);

    virtual ~StateSpaceGrid();

    virtual index_t size() const override;
    virtual const CoordinateSystem& crs() const override;

    virtual const AABox& extent() const override;
    virtual const ArrayShape& shape() const override;
    //virtual const Array& mask() const = 0;

private:
    struct Data;
    std::unique_ptr<Data> mData;

};




/**
 * State space partitioning scheme that operates on gridded state spaces.
 * 
 * The partitioning scheme divides the state space into rectangular local domains of fixed size.
 * Domains may be fully disjoint or partly overlap. In the latter case, the overlapping regions of
 * adjacent local domains are blended together to remove any visible boundaries in the global
 * state after observations have been assimilated.
 * 
 */ 
class GridStateSpacePartitioning : public StateSpacePartitioning
{
public:

    GridStateSpacePartitioning(std::shared_ptr<const GriddedStateSpace> stateSpace, 
                               int blockSize, int padding = 0); 
                            
    virtual ~GridStateSpacePartitioning();

    virtual int numDomains() const override;
    virtual index_t getLocalStateSize(int d) const override;
    virtual void getLocalState(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const override;
    virtual void putLocalState(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;

};








}

#endif