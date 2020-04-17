/**
 * @file Random.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_RANDOM_RANDOM_HPP__
#define __ENDAS_RANDOM_RANDOM_HPP__

#include <Endas/Config.h>
#include <random>

namespace endas
{

/** Random number generator engine used by EnDAS. */
typedef std::mt19937 rng_engine_t;


/**
 * Returns the random number generation engine used by EnDAS.
 */
ENDAS_DLL rng_engine_t& getRngEngine();


}

#endif