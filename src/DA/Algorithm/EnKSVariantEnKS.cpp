#include <Endas/DA/Algorithm/EnsembleKalmanSmoother.hpp>
#include <Endas/Core/Ensemble.hpp>
#include <Endas/Endas.hpp>

#include "../../Compatibility.hpp"

using namespace std;
using namespace endas;


unique_ptr<EnKSVariant> EnKS::clone() const 
{ 
    return make_unique<EnKS>(*this); 
}


void EnKS::processGlobalEnsemble(const Ref<const Array2d> Eg, const ObservationOperator& H,
                                 int k, vector<Array2d>& dataOut) const
{
    dataOut.reserve(2);

    int n = Eg.rows();
    int N = Eg.cols();
    int nobs = H.nobs();

    Array xg = ensembleMean(Eg); 

    // H(Eg - xg) -> dataOut[0]
    dataOut.emplace_back(nobs, N);

    H.apply(Eg.colwise() - xg, dataOut[0]);
    
    // H(Eg) -> dataOut[1]
    dataOut.emplace_back(nobs, N);
    H.apply(Eg, dataOut[1]);
    
}

void EnKS::ensembleTransform(Ref<Array2d> E, vector<Array2d>& Egdata, 
                             const Ref<const Array> z, const CovarianceOperator& R, int k, 
                             Ref<Matrix> Xout) const
{
    ENDAS_ASSERT(Egdata.size() == 2);
    auto HX = Egdata[0].matrix();  
    auto HE = Egdata[1].matrix();

    int n = E.rows();
    int N = E.cols();
    int nobs = HE.rows();

    Matrix K;

    // We have many observations or the observation error covariance operator does not support direct 
    // addition -> perform inversion using the ensemble approximation to R
    if (false && R.mcOnly())
    {
        ENDAS_NOT_IMPLEMENTED;
    }
    // R can be included explicitly in the inversion
    else
    {
        Matrix F = HX * HX.transpose();
        R.fmadd(F, N-1);  // F = F + (N-1)R

        ENDAS_PERF_BEGIN(Invert);
        Eigen::LLT<Ref<Matrix>> cholF(F);
        K = cholF.solve(HX.matrix()).transpose();     // N x nobs
        ENDAS_PERF_END(Invert);
    }

    
    // Sample random perturbations and add them to z
    auto& D = HX;  // We can reuse HX because we no longer need it and it is of correct size
    R.randomMultivariateNormal(D);    
    toAnomaly(D, D); // Ensure zero mean!
    D.colwise()+= z.matrix();

    D-= HE;  // At this point D = (z + N(0, R)) - H(E)

    Xout.noalias() = K * D;
    Xout.diagonal(0).array() += 1.0;

    // Have X transform, apply it to the ensemble
    E.matrix() = E.matrix() * Xout;


}

