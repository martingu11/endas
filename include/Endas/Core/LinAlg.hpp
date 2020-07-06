/**********************************************************************
 * EnDAS - Ensemble Data ASsimilation library
 * Copyright (c) 2020 EnDAS developers
 * 
 * EnDAS is released under the MIT license, see LICENSE.txt.
 **********************************************************************
 * 
 * @file LinAlg.hpp
 * @author Martin Gunia
 * 
 * Basic linear algebra types and operations not offered by Eigen.
 * 
 **********************************************************************/

#ifndef __ENDAS_LINALG_HPP__
#define __ENDAS_LINALG_HPP__

#include <stddef.h>
#include <vector>

#include <Endas/Config.h>
#include <Eigen/Dense>

namespace endas
{

/** 
 * @addtogroup core
 * @{ 
 */


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


/** 
 * An array of indexes.
 * This is aimed at subsetting arrays and matrices. 
 * 
 * @rst
 * .. note::
 *    IndexArray is a dynamic array and is an alias to `std::vector` rather than `Eigen::Array`. 
 *    This is because index arrays are most commonly constructed without prior knowledge of the 
 *    size.
 * @endrst
 */
typedef std::vector<index_t> IndexArray;


/**
 * One-dimensional array of complex numbers. 
 */
typedef Eigen::Array<std::complex<real_t>, Eigen::Dynamic, 1> ComplexArray;

/**
 * Two-dimensional array of complex numbers. 
 */
typedef Eigen::Array<std::complex<real_t>, Eigen::Dynamic, Eigen::Dynamic> ComplexArray2d;



/** Matrix type. */
typedef Eigen::MatrixXd Matrix;

/* Column vector type. Not to be used in the public API, mostly for internal use. */
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1> ColVec;


/** 
 * Opaque reference to a matrix or array object (actual matrix, array or expressions) with 
 * contiguous memory layout.
 * 
 * @rst
 * .. note::
 *    Ref is tailored for performance as most arrays passed to functions are expected to be 
 *    **contiguous**. This enables Eigen to vectorize some expressions, avoids additional pointer
 *    arithmetic on element access and prevents cache trashing. The downside is that when 
 *    non-contiguous arrays are passed to a const Ref, temporary (contiguous) **copy is made**. 
 *    Passing non-contiguous arrays via non-const Ref will result in compilation failure. To avoid
 *    both scenarios, use SoftRef where needed.
 * @endrst
 * 
 */
template <class Derived> using Ref = Eigen::Ref<Derived>;

/**
 * Opaque reference to a matrix or array object (actual matrix, array or expressions).
 * Unlike Ref, SoftRef does not require the referenced memory layout to be contiguous along either 
 * dimension.
 * 
 * This should be used only where actually needed, i.e. when it is expected that non-contiguous
 * arrays will be passed. See Ref documentation for more info.
 */
template <class Derived> using SoftRef = Eigen::Ref<Derived, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;


/**
 * Shape of a multi-dimensional array.
 */
typedef Eigen::Array<index_t, Eigen::Dynamic, 1> ArrayShape;

/**
 * Shape of a two-dimensional array.
 */
typedef Eigen::Array<index_t, 2, 1> ArrayShape2d;

/**
 * Shape of a three-dimensional array.
 */
typedef Eigen::Array<index_t, 3, 1> ArrayShape3d;


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



ENDAS_DLL ArrayShape makeShape(index_t a);
ENDAS_DLL ArrayShape makeShape(index_t a, index_t b);
ENDAS_DLL ArrayShape makeShape(index_t a, index_t b, index_t c);
ENDAS_DLL ArrayShape makeShape(index_t a, index_t b, index_t c, index_t d);
ENDAS_DLL ArrayShape2d makeShape2d(index_t a, index_t b);
ENDAS_DLL ArrayShape3d makeShape3d(index_t a, index_t b, index_t c);




/**
 * Selects array elements `A[i]` for all `i` in `indices`.
 * The selected subset is written into the `out` array which should be pre-allocated to the 
 * correct size (`indices.size()`). 
 * 
 * @param A         The input array to select elements from
 * @param indices   Array of element indices 
 * @param out       Array where the selected elements are copied
 */ 
ENDAS_DLL void select(const Ref<const Array> A, const IndexArray& indices, Ref<Array> out);

/**
 * Selects rows `A[i,:]` for all `i` in `rows`.
 * The selected subset is written into the `out` array which should be pre-allocated to the 
 * correct size (`rows.size()` x `A.cols()`). 
 * 
 * @param A     The input array to select rows from
 * @param rows  Array of row indices 
 * @param out   Array where the selected rows are copied
 */ 
ENDAS_DLL void selectRows(const Ref<const Array2d> A, const IndexArray& rows, Ref<Array2d> out);


/**
 * Selects columns `A[:,j]` for all `j` in `cols`.
 * The selected subset is written into the `out` array which should be pre-allocated to the 
 * correct size (`A.rows()` x `cols.size()`). 
 *  
 * @param A     The input array to select columns from
 * @param cols  Array of column indices 
 * @param out   Array where the selected columns are copied
 */ 
ENDAS_DLL void selectCols(const Ref<const Array2d> A, const IndexArray& cols, Ref<Array2d> out);


/**
 * Selects rows and columns `A[i,j]` for all `i` in `rows` and `j` in `cols`.
 * The selected subset is written into the `out` array which should be pre-allocated to the 
 * correct size (`rows.size()` x `cols.size()`).
 * 
 * @param A     The input array to select rows and columns from
 * @param rows  Array of row indices 
 * @param cols  Array of column indices 
 * @param out   Array where the selected rows and columns are copied  
 */ 
ENDAS_DLL void selectRowsCols(const Ref<const Array2d> A, const IndexArray& rows, 
                              const IndexArray& cols, Ref<Array2d> out);

/** 
 * Copies rows of the array `A` to corresponding rows in `out`.
 * This is the opposite operation to selectRows().
 * 
 * @param A     The input array to select rows and columns from
 * @param rows  Array of row indices 
 * @param cols  Array of column indices 
 * @param out   Array where the selected rows and columns are copied  
 */ 
ENDAS_DLL void distributeRows(const Ref<const Array2d> A, const IndexArray& rows, Ref<Array2d> out);

/** 
 * Copies columns of the array `A` to corresponding columns in `out`.
 * This is the opposite operation to selectRows().
 */ 
ENDAS_DLL void distributeCols(const Ref<const Array2d> A, const IndexArray& cols, Ref<Array2d> out);


/*inline Eigen::Map<Array> reshaped(const ArrayBase<Array2d> A)
{
#if EIGEN_HAS_RESHAPE
    return 
#else
    return Eigen::Map<Array
#endif
}*/



/**
 * Calls `fn(x)` for each element `x` in `A`. 
 */
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


/**
 * Calls `fn(x)` for each element `x` in `A`. 
 */
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
 * Reads array data from file stored in NumPy binary format (`.npy`).
 * 
 * Please note that the functionality is currently rather limited. Only 1 and 2-dimensional 
 * double-precision real arrays are supported. If the NumPy array is in row-major (C) order, 
 * it is transposed automatically.
 */
ENDAS_DLL Array2d loadFromNpy(std::string path);

/**
 * Saves the contents of the array to disk in NumPy binary format (`.npy`).
 */
ENDAS_DLL void saveAsNpy(const Ref<const Array2d> A, std::string path);



/**
 * Computes inverse symmetric square root of A. 
 */
ENDAS_DLL void inverseSymmetricSqrt(const Ref<const Matrix> A, Ref<Matrix> out, bool noalias = false);




/** @} */

}

#endif