#include <Endas/Parallel/SerialJobExecutor.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"

#include <ctpl_stl.h>
#include <deque>

using namespace std;
using namespace endas;


struct SerialJobExecutor::Data
{
    std::deque<unique_ptr<AsyncJob>> scheduled;
};


SerialJobExecutor::SerialJobExecutor()
: mData(make_unique<Data>())
{ }

SerialJobExecutor::~SerialJobExecutor() 
{ }

void SerialJobExecutor::beginJobs(int flags) const 
{
    // Flags are ignored since there is nothing to do for singlethreaded execution
}

int SerialJobExecutor::maxConcurrency() const { return 1; }


void SerialJobExecutor::enqueue(unique_ptr<AsyncJob> job) const
{
    job->run(0);
    mData->scheduled.push_back(move(job));

}

void SerialJobExecutor::waitAllCompleted(function<void(AsyncJob&)> onCompleted) const
{
    while (mData->scheduled.size() > 0)
    {
        unique_ptr<AsyncJob> job = move(mData->scheduled.front());
        mData->scheduled.pop_front();

        ENDAS_ASSERT(job);
        onCompleted(*job);
    }
}



void SerialJobExecutor::pipeline(function<unique_ptr<AsyncJob>(int)> source,
                                 function<void(AsyncJob&)> sink,
                                 int maxJobs, int flags) const
{
    unique_ptr<AsyncJob> job;
    int id = 0;

    while ((job = source(id++)) != nullptr)
    {
        job->run(0); // 'Thread ID' will be always 0
        if (sink) sink(*job);
    }
}
