/**
 * @file MemoryArrayCache.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_MEMORY_ARRAY_CACHE_HPP__
#define __ENDAS_MEMORY_ARRAY_CACHE_HPP__

#include "ArrayCache.hpp"


namespace endas
{

/**
 * Trivial array cache implementation relying entirely on main memory.
 */
class ENDAS_DLL MemoryArrayCache : public ArrayCache
{
public:

    MemoryArrayCache();
    MemoryArrayCache(const MemoryArrayCache&) = delete;
    virtual ~MemoryArrayCache();

    virtual handle_t put(const Ref<const Array2d> data) override;
    virtual std::shared_ptr<Array2d> get(handle_t handle) const override;
    virtual void remove(handle_t handle) override;
    virtual void clear() override;
    virtual void markDirty(handle_t handle);

private:
    struct Data;
    std::unique_ptr<Data> mData;
};




}

#endif