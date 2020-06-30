/**
 * @file ArrayIO.hpp
 * @author Martin Gunia
 * 
 */
#ifndef __ENDAS_IO_ARRAYIO_HPP__
#define __ENDAS_IO_ARRAYIO_HPP__

#include <Endas/Core/LinAlg.hpp>

namespace endas
{

/** 
 * @addtogroup io
 * @{ 
 */

/**
 * Reads array data from file stored in NumPy binary format (`.npy`).
 * 
 * Please note that the functionality is currently rather limited. Only 1 and 2-dimensional 
 * double-precision real arrays are supported. If the NumPy array is in row-major (C) order, 
 * it is transposed automatically.
 */
ENDAS_DLL Array2d loadArrayFromNpy(std::string path);

/**
 * Saves the contents of the array to disk in NumPy binary format (`.npy`).
 */
ENDAS_DLL void saveArrayAsNpy(const Ref<const Array2d> A, std::string path);



/** @} */

}

#endif