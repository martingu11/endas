/**
 * @file ArrayCache.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_ARRAY_CACHE_HPP__
#define __ENDAS_ARRAY_CACHE_HPP__

#include <Endas/Core/LinAlg.hpp>
#include <memory>


namespace endas
{

/**
 * Abstract array data cache.
 */
class ENDAS_DLL ArrayCache
{
public:


    /** Handle type used for representing stored arrays. */
    typedef int handle_t;

    static const handle_t NullHandle = -1;


    virtual ~ArrayCache();
    /** 
     * Places an object into the cache and returns handle for retrieval. The data is always 
     * copied.
     * 
     * @param data      Array or matrix to be stored.
     * 
     * @return Handle to the object that represents the stored array.
     */
    virtual handle_t put(const Ref<const Array2d> data) = 0;


    /**
     * Retrieves data from the cache with exclusive access.
     * 
     * @param handle    Handle to the array instance to retrieve.
     * @return  Shared pointer holding the array.
     * 
     * The caller is granted exclusive access to the returned array during the time he/she
     * holds the shared pointer. During that time the array cannot be reclaimed (for example
     * swapped to disk to make space) by the cache. If the array has been modified, markDirty()
     * must be called to notify the cache that the array contents have changed (which may
     * require new sync to disk, for example).
     */
    virtual std::shared_ptr<Array2d> get(handle_t handle) const = 0;


    /**
     * Retrieves data with exclusive access and removes it from the cache.
     *
     * This is a convenience operation which is equaivalent to calling get() followed by
     * remove(). The returned array will be removed form the cache but will only be 
     * destroyed when there are no references to it.
     * 
     * @param handle    Handle to the array instance to retrieve.
     *
     * @return  Original data object.
     */
    virtual std::shared_ptr<Array2d> pop(handle_t handle);


    /** 
     * Removes array from the cache.
     * 
     * @param handle    Handle to the array instance to remove.
     */
    virtual void remove(handle_t handle) = 0;

    /**
     * Removes all data form the cache. This also invalidates all previously returned handles.
     */
    virtual void clear() = 0;

    /**
     * Marks array as modified.
     * The array should be checked out via get(), otherwise the operation does nothing.
     */
    virtual void markDirty(handle_t handle) = 0;

};




}

#endif