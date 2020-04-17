#ifndef __ENDAS_OBSERVATION_OPERATOR_HPP__
#define __ENDAS_OBSERVATION_OPERATOR_HPP__

#include <Endas/Core/LinAlg.hpp>

namespace endas
{

class ENDAS_DLL ObservationOperator
{
public:

    virtual ~ObservationOperator();

    virtual size_t nstate() const = 0;
    virtual size_t nobs() const = 0;

    virtual bool isLinear() const = 0; 

    virtual bool toMartix(Matrix out) const = 0;

};




}

#endif