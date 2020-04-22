#include <Endas/DA/StateSpaceGeneric.hpp>
#include "../Compatibility.hpp"


using namespace std;
using namespace endas;



Generic1dStateSpace::Generic1dStateSpace(index_t size)
: mSize(size), mCRS(1)
{ }

index_t Generic1dStateSpace::size() const 
{
    return mSize;
}

const CoordinateSystem& Generic1dStateSpace::crs() const
{
    return mCRS;
}

 
GenericStateSpacePartitioning::GenericStateSpacePartitioning(std::shared_ptr<const StateSpace> ss)
: StateSpacePartitioning(ss)
{ }


int GenericStateSpacePartitioning::numDomains() const { return mSS->size(); };


index_t GenericStateSpacePartitioning::getLocalStateSize(int d) const
{
    return 1;
}

void GenericStateSpacePartitioning::getLocalState(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const
{
    outXl = Xg.row(d);
}

void GenericStateSpacePartitioning::putLocalState(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const
{
    Xg.row(d) = Xl;
}

