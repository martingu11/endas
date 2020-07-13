#include <Endas/Parallel/Parallel.hpp>
#include <Endas/Parallel/ThreadJobExecutor.hpp>
#include <Endas/Parallel/SerialJobExecutor.hpp>

#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;



static shared_ptr<const AsyncJobExecutor> globJobExecutor = nullptr; 



const AsyncJobExecutor& endas::getDefaultJobExecutor()
{
    if (!globJobExecutor)
        //globJobExecutor = make_shared<ThreadJobExecutor>();
        globJobExecutor = make_shared<SerialJobExecutor>();

    return *globJobExecutor;
}


void endas::setDefaultJobExecutor(shared_ptr<const AsyncJobExecutor> executor)
{
    globJobExecutor = executor;
};
