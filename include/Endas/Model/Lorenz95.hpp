/**
 * @file Lorenz95.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_MODELS_LORENZ95_HPP__
#define __ENDAS_MODELS_LORENZ95_HPP__

#include <Endas/Algorithm/Algorithm.hpp>
#include <memory>


namespace endas
{

/**
 * Lorenz 95 evolution model.
 * 
 * The Lorenz 95 model is a dynamical system formulated by Edward Lorenz in [1]. 
 * 
 * @see Lorenz, Edward (1996). "Predictability – A problem partly solved"
 */ 
class ENDAS_DLL Lorenz95Model : public LinearizedEvolutionModel
{
public:

    Lorenz95Model(int n = 40, int F = 8);
    virtual ~Lorenz95Model();

    virtual void apply(Ref<Array2d> x, int k, double dt, bool store = true) const override;
    virtual void tl(Ref<Array2d> x, int k) const override;
    virtual void adj(Ref<Array2d> x, int k) const override;
    virtual void stepFinished(int k) const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};



}


#endif