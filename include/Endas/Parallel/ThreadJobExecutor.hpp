/**
 * @file ThreadJobExecutor.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_THREAD_JOB_EXECUTOR_HPP__
#define __ENDAS_THREAD_JOB_EXECUTOR_HPP__

#include <stddef.h>

#include <Endas/Parallel/AsyncJobExecutor.hpp>


namespace endas
{

/** 
 * @addtogroup parallel
 * @{ 
 */



/**
 * Asynchronous job executor implemented using a thread pool.
 * 
 */
class ENDAS_DLL ThreadJobExecutor : public AsyncJobExecutor
{
public:

    /**
     * Creates new ThreadJobExecutor with default thread pool size.
     * The size of the thread pool will match the number of CPU cores found on the system.
     */
    ThreadJobExecutor();

    /**
     * Creates new ThreadJobExecutor with given thread pool size.
     */
    ThreadJobExecutor(int maxThreads);

    ~ThreadJobExecutor();

    /**
     * Resizes the internal thread pool.
     * 
     * The method can only be called when no jobs are currently scheduled.
     */
    void setMaxThreads(int maxThreads);

    virtual int maxConcurrency() const override;
    virtual void beginJobs(int flags = AEFDefault) const override;
    virtual void enqueue(std::unique_ptr<AsyncJob> job) const override;
    virtual void waitAllCompleted(std::function<void(AsyncJob&)> onCompleted) const override;

    virtual void pipeline(std::function<std::unique_ptr<AsyncJob>(int)> source,
                          std::function<void(AsyncJob&)> sink,
                          int maxJobs = 0,
                          int flags = AEFDefault) const override;    

private:

    struct Data;
    std::unique_ptr<Data> mData;
};





/** @} */

}

#endif