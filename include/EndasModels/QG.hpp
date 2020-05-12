/**
 * @file QG.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_MODELS_QG_HPP__
#define __ENDAS_MODELS_QG_HPP__

#include <Endas/DA/Model.hpp>
#include <memory>


namespace endas
{

/**
 * 1.5 -layer Quasi-Geosptrophic ocean circulation model.
 * 
 */ 
class ENDAS_DLL QGModel : public EvolutionModel
{
public:

    QGModel(int N);
    ~QGModel();

    void init(const Ref<const Array2d> x0); 


    void calc_psi(const Ref<const Array2d> E, Ref<Array2d> Eout);

    virtual void operator()(Ref<Array2d> x, int k, double dt, bool store = true) const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};



}


#endif