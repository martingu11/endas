
#include <Endas/Caching/ArrayCache.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;

 ArrayCache::~ArrayCache() { }


shared_ptr<Array2d> ArrayCache::pop(handle_t handle)
{
    shared_ptr<Array2d> A = get(handle);
    if (A) remove(handle);
    return A;
}