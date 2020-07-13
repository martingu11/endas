/**
 * @file SerialJobExecutor.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_SERIAL_JOB_EXECUTOR_HPP__
#define __ENDAS_SERIAL_JOB_EXECUTOR_HPP__

#include <stddef.h>

#include <Endas/Parallel/AsyncJobExecutor.hpp>


namespace endas
{

/** 
 * @addtogroup parallel
 * @{ 
 */



/**
 * Trivial job executor that executes jobs serially in the calling thread.
 * 
 * The executor is useful for debugging when paralellism is not desirable or when
 * serial exectution should be exforced for any other reason.
 */
class ENDAS_DLL SerialJobExecutor : public AsyncJobExecutor
{
public:

    SerialJobExecutor();
    ~SerialJobExecutor();

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