#ifndef __ENDAS_ALGORITHMS_KF_HPP__
#define __ENDAS_ALGORITHMS_KF_HPP__


#include <Endas/Model/Model.hpp>
#include <Endas/Algorithm/Algorithm.hpp>
#include <Endas/Caching/ArrayCache.hpp>

#include <memory>

namespace endas
{


/**
 * Full-rank Kalman Filter and Smoother.
 */
class ENDAS_DLL KalmanSmoother : public SequentialSmoother
{
public:

    /** Evolution model required by this algorithm. */
    typedef LinearizedEvolutionModel Model;


    /** 
     * Constructs KalmanSmoother as a sequential filter or smoother.
     * 
     * @param model State evolution model instance.
     * @param lag   Smoother lag (number of time steps).
     * @param cache Array cache instance for storing intermediate filtering solutions. If `nullptr` 
     *              is given, arrays are cached in memory.
     * 
     * If `lag` is 0, only Kalman Filter solutions will be computed. Array caching is not used when 
     * only computing filter solutions. The reference to the model is stored internally and it is up 
     * to the caller to guarantee that the instance remains alive. 
     */
    KalmanSmoother(const Model& model, int lag = 0, std::shared_ptr<ArrayCache> cache = nullptr);

    ~KalmanSmoother();

    virtual void beginSmoother(const Ref<const Array>x0, const Ref<const Matrix> P0, int k0) override;
    virtual void forecast(Ref<Array> x, Ref<Matrix> P, const Ref<const Matrix> Q, int k, double dt) override;
    virtual void beginAnalysis(Ref<Array> x, Ref<Matrix> P, int k) override;
    virtual void assimilate(const Ref<const Array> z, const Ref<const Matrix> H, const Ref<const Matrix> R) override;
    virtual void endAnalysis() override;
    virtual void endSmoother() override;

private:
    struct Data;
    std::unique_ptr<Data> mData;
};



}


#endif
