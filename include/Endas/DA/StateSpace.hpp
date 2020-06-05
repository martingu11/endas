/**
 * @file StateSpace.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_DA_STATESPACE_HPP__
#define __ENDAS_DA_STATESPACE_HPP__

#include <Endas/Core/Core.hpp>
#include <Endas/Core/AABox.hpp>
#include <Endas/Spatial/CoordinateSystem.hpp>

namespace endas
{

/** 
 * @addtogroup da
 * @{ 
 */



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

    /** Rectangular region within the state space grid. */
    typedef Eigen::AlignedBox<index_t, Eigen::Dynamic> Block;


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
     * Returns size of the subset of the state vector that is inside the given rectangular region.
     */
    virtual index_t size(const Block& block) const = 0;


    /**
     * Returns indices to state vector elements within the given block.
     * 
     * @rst
     * .. tip::
     *    Use :func:`hasEfficientSubset()` to determine whether :func:`getIndices()` is more efficient 
     *    than :func:`getSubset()` and :func:`putSubset()`.
     * @endrst
     */
    virtual void getIndices(const Block& block, IndexArray& out) const = 0;


    /** 
     * Returns `true` if the getSubset() method offers more efficient alternative to getIndices().
     * 
     * For dense grids, access to sub-regions via getSubset() and putSubset() is typically more 
     * efficient than subsetting state vector using indices returned from getIndices(). For sparse
     * grids or grids whose state vector is shuffled in some way the getIndices() approach will be
     * more efficient. 
     * 
     * The default implementation returns ``false``.
     */
    virtual bool hasEfficientSubset() const;

    /** 
     * Reads subset of the state vector for given block.
     * 
     * @rst
     * .. tip::
     *    The default implementation falls back on calling :func:`getIndices()` and will perform poorly if 
     *    :func:`putSubset()` is also called at some point (because indices are computed twice). Use 
     *    :func:`hasEfficientSubset()` to determine whether :func:`getSubset()` and :func:`putSubset()` 
     *    offer efficient access and use :func:`getIndices()` directly otherwise (likely storing the 
     *    indices for later as well).
     * @endrst
     * 
     * @param block  Subset of the grid to read.
     * @param X      State vector or an ensemble of state vectors to read from.
     * @param out    State vector or an ensemble of state vectors to write the subset to.
     */
    virtual void getSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const;


    /** 
     * Writes subset of the state vector for given block.
     *
     * See getSubset() documentation for more information on performance of getSubset() and putSubset().
     * 
     * @param block  Subset of the grid to write.
     * @param X      State vector or an ensemble of state vectors to read from.
     * @param out    State vector or an ensemble of state vectors to write to.
     */
    virtual void putSubset(const Block& block, const Ref<const Array2d> X, Ref<Array2d> out) const;


    /**
     * Returns flat array of indexes identifying grid cells included in the state vector. 
     * 
     * If empty array is returned (`size()==0`), all grid cells are assumed to be included
     * in the state vector.
     */
    //virtual const Array& mask() const = 0;

};


/** @} */

}


#endif