/**
 * @file AABox.hpp
 * @author Martin Gunia
 * 
 * Axis-aligned box.
 */

#ifndef __ENDAS_CORE_AABOX_HPP__
#define __ENDAS_CORE_AABOX_HPP__

#include <Endas/Config.h>
#include <Eigen/Geometry>


namespace endas
{

/** 
 * @addtogroup core
 * @{ 
 */


/**
 * Axis-aligned box with dynamic dimension.
 */
typedef Eigen::AlignedBox<double, Eigen::Dynamic> AABox;


/**
 * Axis-aligned box in one dimension.
 */
typedef Eigen::AlignedBox<double, 1> AABox1d;

/**
 * Axis-aligned box in two dimensions.
 */
typedef Eigen::AlignedBox<double, 2> AABox2d;

/**
 * Axis-aligned box in three dimensions.
 */
typedef Eigen::AlignedBox<double, 3> AABox3d;


/** Rectangular region in a discrete space N-dimensional space. */
typedef Eigen::AlignedBox<index_t, Eigen::Dynamic> IntBox;

/** Rectangular region in a discrete space 1-dimensional space. */
typedef Eigen::AlignedBox<index_t, 1> IntBox1d;

/** Rectangular region in a discrete space 2-dimensional space. */
typedef Eigen::AlignedBox<index_t, 2> IntBox2d;

/** Rectangular region in a discrete space 3-dimensional space. */
typedef Eigen::AlignedBox<index_t, 3> IntBox3d;


template <class Box> Box makeBox(typename Box::Scalar a, typename Box::Scalar b)
{
    typename Box::VectorType _min(1);
    typename Box::VectorType _max(1);
    _min << a;
    _max << b;
    return Box(_min, _max);
}

template <class Box> Box makeBox(typename Box::Scalar a1, typename Box::Scalar a2,
                                 typename Box::Scalar b1, typename Box::Scalar b2)
{
    typename Box::VectorType _min(2);
    typename Box::VectorType _max(2);
    _min << a1, a2;
    _max << b1, b2;
    return Box(_min, _max);
}

template <class Box> Box makeBox(typename Box::Scalar a1, typename Box::Scalar a2, typename Box::Scalar a3,
                                 typename Box::Scalar b1, typename Box::Scalar b2, typename Box::Scalar b3)
{
    typename Box::VectorType _min(3);
    typename Box::VectorType _max(3);
    _min << a1, a2, a3;
    _max << b1, b2, b3;
    return Box(_min, _max);
}



/** @} */

}

#endif