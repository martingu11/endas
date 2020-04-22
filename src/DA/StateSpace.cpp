#include <Endas/DA/StateSpace.hpp>

using namespace std;
using namespace endas;



StateSpace::~StateSpace()
{ }

int GriddedStateSpace::ndim() const
{
    return this->crs().ndim();
}


StateSpacePartitioning::StateSpacePartitioning(shared_ptr<const StateSpace> ss)
: mSS(ss)
{ 
    ENDAS_ASSERT(ss);
}

const StateSpace& StateSpacePartitioning::stateSpace() const
{
    return *mSS;
}

StateSpacePartitioning::~StateSpacePartitioning()
{ }