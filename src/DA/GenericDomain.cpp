#include <Endas/DA/GenericDomain.hpp>
#include "../Compatibility.hpp"


using namespace std;
using namespace endas;


// PartitionPointQuery that works with GenericDomain.
// GenericDomain is a one-dimensional space with each element's index being also its coordinate.
class GenericSSPointQuery : public PartitionPointQuery
{
public:

    GenericSSPointQuery(int stateSize, Array2d coords) 
    : mSize(stateSize), mIndexedCoords(move(coords))
    { 
        ENDAS_ASSERT(coords.rows() == 1 && 
            "Observation coordinate array does not match the dimension of the domain partitioner");
    }

    virtual void rangeQuery(int domain, double range, IndexArray& out) const
    {
        // The domain ID = state variable = point coordinate so the search is trivial.
        // This is meant for small problems so exhaustive search is just fine

        int start = domain - (int)range;
        int end = domain + (int)range;

        for (int i = 0; i != mIndexedCoords.size(); i++) 
        {
            int coord = (int)mIndexedCoords(i);
            if (coord >= start && coord <= end) out.push_back(i);
        }
    }

private:
    int mSize;
    Array2d mIndexedCoords;  
};



GenericDomain::GenericDomain(index_t size)
: mSize(size)
{ }

index_t GenericDomain::size() const { return mSize; }
int GenericDomain::partitionCoordDim() const { return 1; }

const DiscreteDomain& GenericDomain::domain() const 
{ 
    return *this; 
}

int GenericDomain::numLocalDomains() const 
{ 
    return mSize; 
};

index_t GenericDomain::getLocalSize(int d) const
{
    return 1;
}

void GenericDomain::getLocal(int d, const Ref<const Array2d> Xg, Ref<Array2d> outXl) const
{
    outXl = Xg.row(d);
}

void GenericDomain::putLocal(int d, const Ref<const Array2d> Xl, Ref<Array2d> Xg) const
{
    Xg.row(d) = Xl;
}


shared_ptr<const PartitionPointQuery> GenericDomain::indexPoints(Array2d coords) const
{
    return make_shared<GenericSSPointQuery>((int)mSize, move(coords));
}
