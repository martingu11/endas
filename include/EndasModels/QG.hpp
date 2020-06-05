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
 * 1.5 -layer Quasi-Geosptrophic (QG) ocean circulation model.
 * 
 * For more details see 
 * 
 * SAKOV, P. and OKE, P.R. (2008), A deterministic formulation of the ensemble Kalman filter: 
 * an alternative to ensemble square root filters. Tellus A, 60: 361-371. 
 * doi:10.1111/j.1600-0870.2007.00299.x
 */ 
class ENDAS_DLL QGModel : public EvolutionModel
{
public:

    /**
     * Creates new QG model instance. 
     * 
     * @param N             Ensemble size 
     * @param inernalStep   Maximum model integration time step used internally by the RK4 scheme
     */
    QGModel(int N, double internalStep = 1e30);

    ~QGModel();

    /** Returns size of the model discretization grid on the x axis. */
    int sizex() const;

    /** Returns size of the model discretization grid on the y axis. */
    int sizey() const;

    virtual void operator()(Ref<Array2d> x, int k, double dt, bool store = true) const override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};



}


#endif