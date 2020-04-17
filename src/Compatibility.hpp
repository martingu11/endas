#ifndef __ENDAS_COMPATIBILITY_H__
#define __ENDAS_COMPATIBILITY_H__

#include <Endas/Config.h>
#include <memory>

namespace endas
{

#if (!ENDAS_COMPILER_IS_MSVC) && (__cplusplus < 201402L)

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

#else 

using std::make_unique;

#endif


}

#endif

