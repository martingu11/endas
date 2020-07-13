/**
 * @file AsyncJobExecutor.hpp
 * @author Martin Gunia
 */
#ifndef __ENDAS_ASYNC_JOB_EXECUTOR_HPP__
#define __ENDAS_ASYNC_JOB_EXECUTOR_HPP__

#include <Endas/Config.h>
#include <functional>

namespace endas
{

/** 
 * @addtogroup parallel
 * @{ 
 */


/**
 * Piece of work executed by AsyncJobExecutor.
 *
 */
class ENDAS_DLL AsyncJob
{
public:

    virtual ~AsyncJob();

    /**
     * Executes the job.
     * 
     * @param id    Unique identifier of the resource executing the job (thread id, ...).
     */
    virtual void run(int id) = 0;

};


/**
 * Flags controlling execution of asynchronous jobs.
 */ 
enum AsyncExecutionFlags
{
    /** No special flags. */
    AEFNone             = 0x00000000,

    /** 
     * Set up Eigen to run single-threaded if that is expected to increase performance.
     * For executors that spawn their own threads, jobs using Eigen operations are better
     * run single-threaded to avoid oversubscription.
     */
    AEFSetEigenThreads  = 0x00000001,


    /** Default execution flags. */
    AEFDefault          = AEFSetEigenThreads
};



/**
 * Asynchronous job executor.
 */
class ENDAS_DLL AsyncJobExecutor
{
    public:

    virtual ~AsyncJobExecutor();


    /**
     * Returns the maximum number of jobs running concurrently.
     * 
     * @rst
     * .. note::
     *    Calling this method after beginJobs() may give more accurate value on some 
     *    implementations. 
     * @endrst
     */
    virtual int maxConcurrency() const = 0;


    /**
     * Signals to the executor that jobs will be enqueued.
     * 
     * This should be called before any jobs are submitted via enqueue() to initialize necessary 
     * resources (threads, MPI communicators). AsyncJobExecutor implementations are encouraged 
     * to postpone creation of these until they are actually needed, that is when a call to 
     * beginJobs() is made.
     */
    virtual void beginJobs(int flags = AEFDefault) const = 0;

    /** 
     * Schedules a single job for execution.
     * 
     * The  ownership of the job is transferred to the executor.
     */
    virtual void enqueue(std::unique_ptr<AsyncJob> job) const = 0;


    /**
     * Waits until all scheduled jobs have completed. 
     * 
     * The `onCompleted` callback is called for every job that has been completed. 
     */
    virtual void waitAllCompleted(std::function<void(AsyncJob&)> onCompleted) const = 0;


    /** 
     * Implements job pipeline. 
     * 
     * The `source` callable is called to produce jobs to be executed in parallel. Once completed, 
     * jobs are passed to the `sink` callback. Both `source` and `sink` callables run on the calling
     * threads while the individual jobs are processed in parallel. Processing ends when `source` 
     * returns `nullptr`. The main difference compared to sumbitting all jobs via enqeue() is that 
     * the number of jobs currently in the pipeline can be explicitly set. Under the default setting,
     * new jobs are only created once there are available resources for executing them. The maximum
     * can however be lowered so that at most `maxJobs` job instances exist simultaneously. This can 
     * be useful when jobs acquire expensive resources and the number of jobs existing at the same 
     * time should therefore be controlled.
     * 
     * @param source    Callable producing jobs to be processed. The only argument is a counter 
     *                  increased after each call, starting at 0
     * @param sink      Callable executed once a job has been completed
     * @param maxJobs   Maximum number of jobs that are allowed to exist simultaneously. Passing 
     *                  zero or negative number instructs the implementation to select a value based
     *                  on available execution resources
     */
    virtual void pipeline(std::function<std::unique_ptr<AsyncJob>(int)> source,
                          std::function<void(AsyncJob&)> sink,
                          int maxJobs = 0, 
                          int flags = AEFDefault) const = 0;




};


/** @} */

}

#endif