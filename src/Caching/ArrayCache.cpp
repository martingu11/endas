
#include <Endas/Caching/ArrayCache.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;



ArrayCacheEntry::ArrayCacheEntry() : isDirty(false) { }

ArrayCacheEntry::~ArrayCacheEntry() { }


ArrayCache::~ArrayCache() { }

shared_ptr<ArrayCacheEntry> ArrayCache::pop(handle_t handle)
{
    shared_ptr<ArrayCacheEntry> A = get(handle);
    if (A) remove(handle);
    return A;
}