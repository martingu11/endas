/**
 * @file LinAlg.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_LINALG_HPP__
#define __ENDAS_LINALG_HPP__

#include <stddef.h>

#include <Endas/Config.h>
#include <Eigen/Dense>

namespace endas
{

/**
 * One-dimensional array type. When used in linear algebra operations, the array is interpreted
 * as a column vector.
 */
typedef Eigen::Array<real_t, Eigen::Dynamic, 1> Array;

/** Two-dimensional array. */
typedef Eigen::Array<real_t, Eigen::Dynamic, Eigen::Dynamic> Array2d;


/** Matrix type. */
typedef Eigen::MatrixXd Matrix;

/* Column vector type. Not to be used in the public API, mostly for internal use. */
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1> ColVec;


/** 
 * Opaque reference to a matrix or array object (actual matrix, array or expressions).
 */
template <class Derived> using Ref = Eigen::Ref<Derived>;

/**
 * Opaque reference to a matrix or array object (actual matrix, array or expressions).
 * Unlike Ref<>, SoftRef<> does not require the referenced memory layout to be contiguous.
 */
template <class Derived> using SoftRef = Eigen::Ref<Derived, 0, Eigen::InnerStride<>>;


/** 
 * Returns reference to a shared instance of an empty matrix object. 
 */
ENDAS_DLL const Matrix& emptyMatrix();


/** 
 * Returns new Array instance populated with values from the initializer list.
 * 
 *     Array a = makeArray({1, 2, 3});
 */
ENDAS_DLL Array makeArray(std::initializer_list<real_t> values);


/** 
 * Returns new Matrix instance populated with values from the initializer list.
 * Regardless of the matrix starage layout, the values should be specified in row-major
 * order in the initializer list (i.e. the way they would be written out).
 * 
 *     Matrix a = makeMatrix(2, 3, {
 *         1, 2, 3,
 *         4, 5, 6
 *     });  
 */
ENDAS_DLL Matrix makeMatrix(int rows, int cols, std::initializer_list<real_t> values);


}

#endif