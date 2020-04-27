
#include <Endas/Random/Random.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"

#include <thread>
#include <functional>
#include <chrono>
#include <mutex>


using namespace std;
using namespace endas;

static shared_ptr<RandomNumberGenerator> globRNG = make_shared<MT19937>();
static mutex globSeedMutex;


RandomNumberGenerator::~RandomNumberGenerator()
{ }


RandomNumberGenerator& endas::getRandomNumberGenerator()
{
    thread_local unique_ptr<RandomNumberGenerator> rng = globRNG->clone();
    ENDAS_ASSERT(rng);
    return *rng;
}

void endas::setRandomNumberGenerator(const RandomNumberGenerator& rng)
{
    globRNG = rng.clone();
}


unsigned int endas::getRandomUniqueSeed()
{
    static random_device rd;
    // Lock here?
    return rd();
}


template <class RNG> StandardRNG<RNG>::StandardRNG()
: mGen((typename RNG::result_type)getRandomUniqueSeed())
{ }


template <class RNG> unique_ptr<RandomNumberGenerator> StandardRNG<RNG>::clone() const
{
    auto clone = make_unique< StandardRNG<RNG> >();
    clone->seed((typename RNG::result_type)getRandomUniqueSeed());
    return move(clone);    
}


template <class RNG> void StandardRNG<RNG>::seed(unsigned int value) 
{
    mGen.seed((typename RNG::result_type)value);
}


template <class RNG> real_t StandardRNG<RNG>::standardNormal()
{
    normal_distribution<real_t> dist;
    return dist(mGen);
}

template <class RNG> void StandardRNG<RNG>::standardNormal(SoftRef<Array2d> out)
{
    normal_distribution<real_t> dist;
    out = out.unaryExpr([&](real_t) { return dist(mGen); });
}



// Note: Declare all StandardRNG instantiations here!

template class StandardRNG<std::mt19937>;


