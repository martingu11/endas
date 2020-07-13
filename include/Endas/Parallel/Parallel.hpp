/**
 * @file Parallel.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_PARALLEL_HPP__
#define __ENDAS_PARALLEL_HPP__

#include <stddef.h>

#include "AsyncJobExecutor.hpp"


namespace endas
{

/** 
 * @addtogroup parallel
 * @{ 
 */


/**
 * Returns the default asynchronous job executor.
 */
ENDAS_DLL const AsyncJobExecutor& getDefaultJobExecutor();


/**
 * Sets the default asynchronous job executor.
 */
ENDAS_DLL void setDefaultJobExecutor(std::shared_ptr<const AsyncJobExecutor> executor);



/** @} */

}

#endif