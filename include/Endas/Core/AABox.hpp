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


/** @} */

}

#endif