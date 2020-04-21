/**
 * @file Model.hpp
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

    template <class CRS>
    StateSpaceGrid(const GridShape& shape, const CRS& crs, const AABox& extent)
    : StateSpaceGrid(shape, std::make_shared<CRS>(crs), extent)
    { }

    StateSpaceGrid(const GridShape& shape, std::shared_ptr<const CoordinateSystem> crs, 
                   const AABox& extent);

    virtual ~StateSpaceGrid();

    virtual index_t size() const override;
    virtual const CoordinateSystem& crs() const override;

    virtual const AABox& extent() const override;
    virtual const GridShape& shape() const override;
    //virtual const Array& mask() const = 0;

private:
    struct Data;
    std::unique_ptr<Data> mData;

};







}

#endif