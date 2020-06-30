/**
 * @file Exception.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_EXCEPTION_HPP__
#define __ENDAS_EXCEPTION_HPP__

#include <Endas/Config.h>
#include <exception>
#include <string>
#include <memory>


namespace endas
{

/** 
 * @addtogroup core
 * @{ 
 */


/**
 * Exception signalling that the requested functionality is not supported.
 */
class ENDAS_DLL NotSupportedError : public std::logic_error
{
public:
    NotSupportedError(const char* why, const char* file, int line);
};

/**
 * Exception signalling that the requested functionality is not yet implemented.
 */
class ENDAS_DLL NotImplementedError : public std::logic_error
{
public:
    NotImplementedError(const char* file, int line);
};

/**
 * Runtime assertion.
 * 
 * EnDAS can be configured to either use `assert()` or throw AssertionError. The latter 
 * will be in use also in release builds and can be used to detect issues when working
 * with large data.
 * 
 * Do not use directly, use the ENDAS_ASSERT macro instead!
 */
class ENDAS_DLL AssertionError : public std::logic_error
{
public:
    AssertionError(const char* expr, const char* file, int line);

};


#define ENDAS_NOT_SUPPORTED(why) { throw NotSupportedError(why, __FILE__, __LINE__); }
#define ENDAS_NOT_IMPLEMENTED { throw NotImplementedError(__FILE__, __LINE__); }


/** 
 * EnDAS assertion.
 */
#ifdef NDEBUG
#   define ENDAS_ASSERT(expr) { if (!(expr)) throw AssertionError(#expr, __FILE__, __LINE__); }
#else 
#   define ENDAS_ASSERT(expr) assert(expr)
#endif 



#define ENDAS_REQUIRE(expr, extype, msg) { if (!(expr)) throw extype(msg); }

#define ENDAS_CHECK_ARGUMENT(expr, msg) { if (!(expr)) throw std::invalid_argument(msg); }



/** @} */

}

#endif