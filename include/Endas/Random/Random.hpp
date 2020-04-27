/**
 * @file Random.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_RANDOM_RANDOM_HPP__
#define __ENDAS_RANDOM_RANDOM_HPP__


#include <Endas/Core/LinAlg.hpp>
#include <random>

namespace endas
{


/**
 * Random number generator (RNG) interface used by EnDAS.
 * 
 * Please note that this is not meant to be a general-purpose random number generator (although it
 * can be used as one if it fits your needs). Instead, RandomNumberGenerator is used by EnDAS to 
 * generate random number sequences while allowing you to supply a RNG of your choice. For this 
 * reason, RandomNumberGenerator only implements sampling from distributions that EnDAS needs.
 */ 
class ENDAS_DLL RandomNumberGenerator
{
public:
    virtual ~RandomNumberGenerator();

    /**
     * Returns new RandomNumberGenerator instance. 
     * 
     * @rst
     * .. attention::
     *    Note to implementers: The clone method must always return random number generator 
     *    initialized from a random seed. Copying the state from ``this`` will result
     *    in identical random number streams generated on different threads, which is almost
     *    definitely not what one wants! You can use :func:`endas::getRandomUniqueSeed()` in 
     *    your implementation to seed the generator.
     * @endrst
     */
    virtual std::unique_ptr<RandomNumberGenerator> clone() const = 0;

    /** 
     * Initializes the internal state of the generator from the given seed. 
     * 
     * @rst
     * .. note::
     *    All generators are initialized with unique seed values by default. Seeding the generator 
     *    explicitly should only be done when testing your DA code or when deterministic output is 
     *    needed for reproducibility. In this case you should also disable multi-threaded execution 
     *    for parts of EnDAS depending on random numbers (via the ``FORCE_DETERMINISTIC`` 
     *    configuration option). 
     * @endrst
     */
    virtual void seed(unsigned int value) = 0;


    /** 
     * Generates a single sample from standard Normal distribution (zero mean and variance 1).
     */
    virtual real_t standardNormal() = 0;

    /** 
     * Fills array with samples drawn from standard Normal distribution (zero mean and variance 1).
     */
    virtual void standardNormal(SoftRef<Array2d> out) = 0;

};


/**
 * Returns the random number generation engine used by EnDAS.
 * 
 * For multi-threaded applications, this returns unique copy of the random number generator for 
 * each calling thread.
 */
ENDAS_DLL RandomNumberGenerator& getRandomNumberGenerator();


/**
 * Installs random number generation engine used by EnDAS.
 * 
 * Please note that the passed instance is cloned, which implies that the internal state of the 
 * passed generator is not carried over. To use a custom generator, setRandomNumberGenerator() 
 * **must** be called before any EnDAS functionality relying on random numbers and likewise before 
 * getRandomNumberGenerator() is called. Any changes made via setRandomNumberGenerator() afterwards
 * will have no effect or may produce inconsitent results. 
 * 
 * Because the internal state does not carry over, you should seed the generator like this::
 * 
 *     endas::setRandomNumberGenerator(your_generator_of_choice);
 *     endas::getRandomNumberGenerator().seed(1234);
 * 
 * @rst
 * .. attention::
 *    All generator are initialized from unique seed values by default. Seeding the generator 
 *    explicitly should only be done when testing your DA code and when deterministic output is 
 *    needed. In this case you must also **disable multi-threaded execution** for EnDAS. This is 
 *    because the generator is internally cloned for each thread and seeded with unique random 
 *    seed (to avoid identical number sequences to be generated in different calling threads).
 * @endrst 
 * 
 */
ENDAS_DLL void setRandomNumberGenerator(const RandomNumberGenerator& rng);


/** 
 * Returns new random seed.
 * The function is thread-safe.
 */
ENDAS_DLL unsigned int getRandomUniqueSeed();


/**
 * RandomNumberGenerator implementation based on C++11 random bit generators.
 * 
 * Currently the following concrete implementations are defined:
 * 
 * - ``endas::MT19937``
 */ 
template <class UniformRandomBitGenerator>
class ENDAS_DLL StandardRNG : public RandomNumberGenerator
{
public:

    StandardRNG();

    virtual std::unique_ptr<RandomNumberGenerator> clone() const override;
    virtual void seed(unsigned int value) override;
    virtual real_t standardNormal() override;
    virtual void standardNormal(SoftRef<Array2d> out) override;

private:
    UniformRandomBitGenerator mGen;
};


/** Random number generator based on the Mersenne Twister algorithm (MT19937 variant). */
typedef StandardRNG<std::mt19937> MT19937;




}

#endif