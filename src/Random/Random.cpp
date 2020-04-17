
#include <Endas/Random/Random.hpp>

using namespace std;
using namespace endas;


rng_engine_t& endas::getRngEngine()
{
    static rng_engine_t gen;
    return gen;
}