#include <Endas/Parallel/ThreadJobExecutor.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"

#include <ctpl_stl.h>
#include <deque>

using namespace std;
using namespace endas;


struct ThreadJobExecutor::Data
{
    int maxThreads;
    ctpl::thread_pool tpool;
    std::deque<std::future<unique_ptr<AsyncJob>>> scheduled;

    int numEigenThreads;
    int execFlags;
};

struct JobWrapper
{
    JobWrapper(unique_ptr<AsyncJob> j) : job(move(j)) { }
    unique_ptr<AsyncJob> job;

    unique_ptr<AsyncJob> operator()(int id) 
    { 
        job->run(id); 
        return move(job);
    }
};



ThreadJobExecutor::ThreadJobExecutor()
: mData(make_unique<Data>())
{
    mData->maxThreads = thread::hardware_concurrency();
}

ThreadJobExecutor::ThreadJobExecutor(int maxThreads)
{
    mData->maxThreads = maxThreads;
}

ThreadJobExecutor::~ThreadJobExecutor()
{ }


void ThreadJobExecutor::setMaxThreads(int maxThreads)
{
    // Only set the flag, we will resize the thread pool on next call to beginJobs()
    mData->maxThreads = maxThreads;
}

int ThreadJobExecutor::maxConcurrency() const { return mData->maxThreads; }


void ThreadJobExecutor::beginJobs(int flags) const 
{
    mData->tpool.resize(mData->maxThreads);

    mData->execFlags = flags;
    mData->numEigenThreads = Eigen::nbThreads();

    // Disable Eigen's threads for jobs
    if (mData->execFlags & AEFSetEigenThreads != 0)
    {
        Eigen::setNbThreads(1);
    }
}


void ThreadJobExecutor::enqueue(unique_ptr<AsyncJob> job) const
{
    ENDAS_ASSERT(mData->tpool.size() > 0 && "Forgot to call beginJobs()?");

    mData->scheduled.push_back(mData->tpool.push(JobWrapper(move(job))));
}


void ThreadJobExecutor::waitAllCompleted(function<void(AsyncJob&)> onCompleted) const
{
    // Check the shceduled jobs completion at the order they were scheduled.
    // Todo: We could use a notify mechanism to check in the order the jobs were actually 
    // completed

    while (mData->scheduled.size() > 0)
    {
        auto& future = mData->scheduled.front();
        ENDAS_ASSERT(future.valid());

        shared_ptr<AsyncJob> job = future.get();
        ENDAS_ASSERT(job);
        onCompleted(*job);

        mData->scheduled.pop_front();
    }

    // Restore Eigen's threads
    if (mData->execFlags & AEFSetEigenThreads != 0)
    {
        Eigen::setNbThreads(mData->numEigenThreads);
    }
}


void ThreadJobExecutor::pipeline(function<unique_ptr<AsyncJob>(int)> source,
                                 function<void(AsyncJob&)> sink,
                                 int maxJobs, int flags) const
{
    this->beginJobs(flags);

    // Todo: This just uses enqueue(). Implement proper pipeline with maxJobs!

    unique_ptr<AsyncJob> job;
    int id = 0;
    while ((job = source(id++)) != nullptr)
    {
        // Schedule the job for execution
        this->enqueue(move(job));
        
    }

    this->waitAllCompleted([&sink](AsyncJob& job)
    {
        if (sink) sink(job);
    });
}
