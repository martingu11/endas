/**
 * @file LinAlg.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_LINALG_HPP__
#define __ENDAS_LINALG_HPP__

#include <stddef.h>
#include <vector>

#include <Endas/Config.h>
#include <Eigen/Dense>

namespace endas
{


/**
 * Basic type for indexing arrays. 
 * This is an alias to `Eigen::Index`.
 */ 
typedef Eigen::Index index_t;


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
 * Shape of a multi-dimensional array.
 */
typedef Eigen::Array<int, Eigen::Dynamic, 1> ArrayShape;

/**
 * Shape of a two-dimensional array.
 */
typedef Eigen::Array<int, 2, 1> ArrayShape2d;

/**
 * Shape of a three-dimensional array.
 */
typedef Eigen::Array<int, 2, 1> ArrayShape3d;



/** 
 * An array of indexes.
 * This is aimed at subsetting arrays and matrices. 
 */
typedef std::vector<index_t> IndexArray;


/** 
 * Returns reference to a shared instance of an empty one-dimensional array object. 
 */
ENDAS_DLL const Array& emptyArray();

/** 
 * Returns reference to a shared instance of an empty two-dimensional array object. 
 */
ENDAS_DLL const Array2d& emptyArray2d();


/** 
 * Returns reference to a shared instance of an empty matrix object. 
 */
ENDAS_DLL const Matrix& emptyMatrix();



/** 
 * Returns new Array instance populated with sequence given by start and end.
 * 
 * @param start     Starting value of the sequence
 * @param end       End value of the sequence (not inclusive)
 * @param step      DIstance between donsecutive sequence elements
 */
ENDAS_DLL Array makeSequence(real_t start, real_t end, real_t step = 1.0);



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



/**
 * Computes matrix product out=AB, where A is an NxM dense matrix and B is an MxM diagonal 
 * matrix. The diagonal matrix B is specified by the array of diagonal elements.
 * 
 * @param A     Dense matrix A
 * @param b     Array of diagonal elements of B
 * @param out   Pre-allocated matrix (same size as `A`) where to store the result. 
 * 
 * @note There are no aliasing issues with this kind of matrix product and `out=A` can be used
 *       to perform the multiplication in-place.
 */
ENDAS_DLL void denseXdiag(const Ref<const Matrix> A, const Ref<const Array> b, Ref<Matrix> out);

/**
 * Computes matrix product out=AB, where A is an NxN diagonal matrix and B is a dense NxM matrix.
 * The diagonal matrix B is specified by the array of diagonal elements.
 * 
 * @param a     Array of diagonal elements of A
 * @param B     Dense matrix B
 * @param out   Pre-allocated matrix (same size as `A`) where to store the result. 
 * 
 * @note There are no aliasing issues with this kind of matrix product and `out=B` can be used
 *       to perform the multiplication in-place.
 */
ENDAS_DLL void diagXdense(const Ref<const Array> a, const Ref<const Matrix> B, Ref<Matrix> out);



ENDAS_DLL void select(const Ref<const Array> A, const IndexArray& indices, Ref<Array> out);

ENDAS_DLL void selectRows(const Ref<const Array2d> A, const IndexArray& rows, Ref<Array2d> out);

ENDAS_DLL void selectCols(const Ref<const Array2d> A, const IndexArray& cols, Ref<Array2d> out);

ENDAS_DLL void selectRowsCols(const Ref<const Array2d> A, const IndexArray& rows, 
                              const IndexArray& cols, Ref<Array2d> out);




/*inline Eigen::Map<Array> reshaped(const ArrayBase<Array2d> A)
{
#if EIGEN_HAS_RESHAPE
    return 
#else
    return Eigen::Map<Array
#endif
}*/




template <class Fn>
inline void foreachCoeff(const Ref<const Array> A, Fn fn)
{
#if EIGEN_HAS_ITERATRORS
    int i = 0;
    for (auto&& x : A) fn(i++, x);
#else
    for (int i = 0; i != A.size(); i++) fn(i, A(i));
#endif
}


template <class Fn>
inline void foreachCoeff(Ref<Array> A, Fn fn)
{
#if EIGEN_HAS_ITERATRORS
    int i = 0;
    for (auto&& x : A) fn(i++, x);
#else
    for (int i = 0; i != A.size(); i++) fn(i, A(i));
#endif
}



/**
 * Computes inverse symmetric square root of A. 
 */
ENDAS_DLL void inverseSymmetricSqrt(const Ref<const Matrix> A, Ref<Matrix> out, bool noalias = false);




}

#endif