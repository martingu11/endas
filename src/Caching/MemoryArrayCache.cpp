
#include <Endas/Caching/MemoryArrayCache.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"

#include <unordered_map>

using namespace std;
using namespace endas;


struct MemoryArrayCache::Data
{
    handle_t handleCounter;
    unordered_map<ArrayCache::handle_t, shared_ptr<Array2d>> entries;
};


MemoryArrayCache::MemoryArrayCache()
: mData(make_unique<Data>())
{ }

MemoryArrayCache::~MemoryArrayCache() 
{ }


ArrayCache::handle_t MemoryArrayCache::put(const Ref<const Array2d> data)
{
    shared_ptr<Array2d> copy = make_shared<Array2d>(data);

    handle_t handle = mData->handleCounter++;
    mData->entries.insert(make_pair(handle, copy));
    return handle;
}


shared_ptr<Array2d> MemoryArrayCache::get(handle_t handle) const
{
    ENDAS_ASSERT(handle != ArrayCache::NullHandle);

    auto it = mData->entries.find(handle);
    ENDAS_ASSERT(it != mData->entries.end());
    return it->second;
}

void MemoryArrayCache::remove(handle_t handle) 
{
    ENDAS_ASSERT(handle != ArrayCache::NullHandle);
    mData->entries.erase(handle);
}

void MemoryArrayCache::clear()
{
    mData->entries.clear();
    mData->handleCounter = 0;
}

void MemoryArrayCache::markDirty(handle_t handle)
{ 
    ENDAS_ASSERT(handle != ArrayCache::NullHandle);
}

