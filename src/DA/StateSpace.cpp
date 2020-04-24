#include <Endas/DA/StateSpacePartitioning.hpp>

using namespace std;
using namespace endas;



StateSpace::~StateSpace()
{ }

int GriddedStateSpace::dim() const
{
    return this->crs().dim();
}

StateSpacePartitioning::~StateSpacePartitioning()
{ }