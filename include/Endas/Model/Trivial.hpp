/**
 * @file Trivial.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_MODELS_TRIVIAL_HPP__
#define __ENDAS_MODELS_TRIVIAL_HPP__

#include "Model.hpp"
#include <memory>


namespace endas
{


/**
 * Trivial state evolution model represented by a matrix.
 */ 
class ENDAS_DLL MatrixModel : public LinearizedEvolutionModel
{
public:

    MatrixModel(const Ref<const Matrix> M);

    MatrixModel(int n, const std::function<void(Ref<Matrix>)> M);

    virtual void apply(Ref<Array2d> x, int k, double dt, bool store = true) const override;
    virtual void tl(Ref<Array2d> x, int k) const override;
    virtual void adj(Ref<Array2d> x, int k) const override;

    const Ref<const Matrix> get() const;

private:
    Matrix mModel;
};





}


#endif