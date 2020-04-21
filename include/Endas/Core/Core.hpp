/**
 * @file Core.hpp
 * @author Martin Gunia
 */

#ifndef __ENDAS_CORE_HPP__
#define __ENDAS_CORE_HPP__

#include <Endas/Config.h>
#include <memory>



namespace endas
{

/** 
 * Exclusive ownership pointer with type erasure.
 * 
 * This is currently just an alias for ` std::shared_ptr<void>` although solution based on
 * unique_ptr with custom deleter would be possible too. But this is low priority.
 */
template <class T> using unique_ptr_ex = std::shared_ptr<T>;


template<typename T, typename ...Args>
unique_ptr_ex<T> make_unique_ex( Args&& ...args )
{
    return std::make_shared<T>( std::forward<Args>(args)... );
}



}

#endif