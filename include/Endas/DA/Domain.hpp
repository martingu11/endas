/**
 * @file Domain.hpp
 * @author Martin Gunia
 * 
 * Abstract space representations.
 */

#ifndef __ENDAS_DA_DOMAIN_HPP__
#define __ENDAS_DA_DOMAIN_HPP__

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
 * Base class for abstract discrete domains.
 * 
 * DiscreteDomain is a base class for representing abstract discretized domains such as 
 * the state or observation space. The only information provided is the size of the space.
 */ 
class DiscreteDomain
{
public:

    virtual ~DiscreteDomain();

    /**
     * Returns the size of the space, i.e. the number of discrete elements making up the space. 
     */
    virtual index_t size() const = 0;
};


/**
 * Abstract representation of discrete domains whose elements can be identified by a spatial 
 * coordinate.
 */ 
class DiscreteSpatialDomain : public DiscreteDomain
{
public:

    /**
     * Returns the number of spatial coordinates used to represent the location of each state
     * vector element. This is a convenient alias for `crs().dim()`.
     */
    int dim() const;

    /**
     * Returns the coordinate system of the state space.
     */
    virtual const CoordinateSystem& crs() const = 0;


    /**
     * Returns spatial coordinates of the elements of the state space. 
     * 
     * The coordinates are stored in the ``out`` array which must be pre-allocated to size 
     * *m* x *n*, where *m* is the number of spatial dimensions as returned by ``dim()`` and *n*
     * is the state vector size.
     */
    virtual void getCoords(Ref<Array2d> out) const = 0;

};




/**
 * Base class for domains with elements organized on a regular (multi-dimensional) grid.
 */
class GriddedDomain : public DiscreteSpatialDomain
{
public:

    /* Rectangular region within the state space grid. */
    typedef IntBox Block;

    /**
     * Returns the spatial extent of the grid.
     */
    virtual const AABox& extent() const = 0;

    /**
     * Returns the shape of the grid (the number of cells in each dimension).
     */
    virtual const ArrayShape& shape() const = 0;


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