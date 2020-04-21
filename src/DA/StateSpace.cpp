#include <Endas/DA/StateSpace.hpp>

using namespace std;
using namespace endas;



StateSpace::~StateSpace()
{ }

int GriddedStateSpace::ndim() const
{
    return this->crs().ndim();
}