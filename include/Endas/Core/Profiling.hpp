/**
 * @file Core.hpp
 * @author Martin Gunia
 * 
 * Simple system for collecting profiling information during algorithm execution.
 */

#ifndef __ENDAS_PROFILING_HPP__
#define __ENDAS_PROFILING_HPP__

#include <Endas/Config.h>

#include <chrono>
#include <ostream>


namespace endas
{

typedef std::chrono::steady_clock perfclock_t;


namespace detail
{

template <class BaseUnit>
inline double elapsedTime(perfclock_t::time_point start, perfclock_t::time_point end)
{
    return std::chrono::duration_cast<BaseUnit>(end - start).count() / 1000.0;
}

#if !ENDAS_PROFILING_DISABLED

ENDAS_DLL void recordTime(const char* key, perfclock_t::time_point start, perfclock_t::time_point end);

struct ENDAS_DLL NewPerfScope
{
    NewPerfScope(const char* key);
    ~NewPerfScope();
};

#endif

}


/** 
 * @addtogroup core
 * @{ 
 */


#define ENDAS_TIMER_BEGIN(timer) endas::perfclock_t::time_point timer_##timer##_begin = endas::perfclock_t::now()
#define ENDAS_TIMER_END(timer) endas::perfclock_t::time_point timer_##timer##_end = endas::perfclock_t::now()

#define ENDAS_TIMER_ELAPSED_SEC(timer) endas::detail::elapsedTime<std::chrono::milliseconds>(timer_##timer##_begin, timer_##timer##_end)
#define ENDAS_TIMER_ELAPSED_MSEC(timer) endas::detail::elapsedTime<std::chrono::microseconds>(timer_##timer##_begin, timer_##timer##_end)


#if !ENDAS_PROFILING_DISABLED

/**
 * Defines new performance measurement scope.
 * The scope exists within the current C++ scope and is automatically removed at the C++ scope exit.
 * Expands to empty operation if `ENDAS_PROFILING_DISABLED` is `true`.
 * 
 * @param name  Name of the scope, must be a valid C++ identifier.
 */
#   define ENDAS_PERF_SCOPE(name) endas::detail::NewPerfScope perfscope_##name(#name)


/**
 * Sets up new performance measurement timer.
 * Expands to empty operation if `ENDAS_PROFILING_DISABLED` is `true`.
 * 
 * @param name  Name of the timer, must be a valid C++ identifier.
 */
#   define ENDAS_PERF_BEGIN(name) endas::perfclock_t::time_point perftimer_##name##_begin = endas::perfclock_t::now()


/**
 * Records the time interval between `ENDAS_PERF_BEGIN(name)` and now and records it into the 
 * current performance measurement scope.
 * Expands to empty operation if `ENDAS_PROFILING_DISABLED` is `true`.
 * 
 * @param name  Name of the timer, must be a valid C++ identifier.
 */
#   define ENDAS_PERF_END(name) endas::detail::recordTime(#name, perftimer_##name##_begin, endas::perfclock_t::now())

#else 
#   define ENDAS_PERF_SCOPE(name) 
#   define ENDAS_PERF_BEGIN(name) 
#   define ENDAS_PERF_END(name) 
#endif




/**
 * Deletes any profiling data collected so far. 
 */
ENDAS_DLL void profilerClear();


/**
 * Prints profiling information to a stream.
 * If profiling is disabled, this prints nothing.
 * 
 * @param os          Output stream to which to print
 * @param maxNesting  Controls how many nested scopes will be printed
 */
ENDAS_DLL void profilingSummary(std::ostream& os, int maxNesting = 5);



/** @} */

}

#endif