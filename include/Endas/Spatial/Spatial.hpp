/**
 * @file Spatial.hpp
 * @author Martin Gunia
 * 
 * Basic spatial types.
 */

#ifndef __ENDAS_SPATIAL_SPATIAL_HPP__
#define __ENDAS_SPATIAL_SPATIAL_HPP__

#include <Endas/Config.h>

#include <Eigen/Geometry>


namespace endas
{

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



}

#endif