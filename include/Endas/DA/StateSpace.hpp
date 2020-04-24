/**
 * @file StateSpace.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_STATESPACE_HPP__
#define __ENDAS_DA_STATESPACE_HPP__

#include <Endas/Core/Core.hpp>
#include <Endas/Spatial/Spatial.hpp>
#include <Endas/Spatial/CoordinateSystem.hpp>


namespace endas
{


/**
 * Base class for state space implementations.
 * 
 * This is a generic state space abstraction that only provides information about the size of the
 * state space.
 * 
 * Please note that StateSpace or derived classes are only needed whe using analysis localization 
 * or with observation managers that spatially map observations to state variables. In other cases 
 * no information needs to be known about the state space.
 */ 
class StateSpace
{
public:

    virtual ~StateSpace();

    /**
     * Returns the size of the state space. 
     */
    virtual index_t size() const = 0;

};



/**
 * Base class for state spaces with elements organized on a regular (multi-dimensional) grid.
 */
class GriddedStateSpace : public StateSpace
{
public:

    /**
     * Returns the grid dimension (1 for discrete line, 2 for rectangular grid etc...).
     * This is a convenient alias for `crs().dim()`.
     */
    int dim() const;

    /**
     * Returns the spatial extent of the grid.
     */
    virtual const AABox& extent() const = 0;

    /**
     * Returns the shape of the grid (the number of cells in each dimension).
     */
    virtual const ArrayShape& shape() const = 0;

    /**
     * Returns the coordinate system of the state space.
     */
    virtual const CoordinateSystem& crs() const = 0;

    /**
     * Returns flat array of indexes identifying grid cells included in the state vector. 
     * 
     * If empty array is returned (`size()==0`), all grid cells are assumed to be included
     * in the state vector.
     */
    //virtual const Array& mask() const = 0;

};




}

#endif