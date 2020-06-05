/**
 * @file Core.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_CORE_HPP__
#define __ENDAS_CORE_HPP__


#include "Exception.hpp"
#include "LinAlg.hpp"
#include <memory>
#include <type_traits>



namespace endas
{

/** 
 * @addtogroup core
 * @{ 
 */

struct NullSharedPtrDeleter
{
    void operator()(const void*) { }
};


/**
 * Wraps reference in a shared pointer for use where shared pointers are required by the API.
 * 
 * @note This simply exposes the reference via a shared pointer with empty deleter. It is the
 * responsibility of the caller to ensure that the referenced object remains alive for as long 
 * as it is used!
 */
template <class T>
inline std::shared_ptr<T> shared_ptr_wrap(T& x)
{
    return std::shared_ptr<T>(&x, NullSharedPtrDeleter());
}


template <class From, class To> 
using require_is_convertible = typename std::enable_if<std::is_convertible<From*, To*>::value, bool>::type;



/** @} */

}

#endif