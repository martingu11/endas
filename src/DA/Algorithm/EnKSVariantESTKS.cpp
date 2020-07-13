#include <Endas/DA/Algorithm/EnsembleKalmanSmoother.hpp>
#include <Endas/Core/Ensemble.hpp>
#include <Endas/Endas.hpp>

#include "../../Compatibility.hpp"

using namespace std;
using namespace endas;


unique_ptr<EnKSVariant> ESTKS::clone() const 
{ 
    return make_unique<ESTKS>(*this); 
}


void ESTKS::init(int n, int N) 
{
    /// Pre-compute the ensemble-space transformation. 
    double a = (1.0 / (double)N) * (1.0 / (1.0 / sqrt(N) + 1.0));

    mT = Matrix::Constant(N, N-1, -a);
    mT.diagonal().array() = 1.0 - a;
    mT.row(N-1).array() = -1.0 / sqrt(N);
}


void ESTKS::applyCovInflation(Ref<Array2d> E, double factor, int k) const
{
    mInflation = factor;
}

void ESTKS::processGlobalEnsemble(const Ref<const Array2d> Eg, const ObservationOperator& H,
                                  int k, vector<Array2d>& dataOut) const
{
    dataOut.reserve(2);

    int n = Eg.rows();
    int N = Eg.cols();

    Array xg = ensembleMean(Eg); 

    // H(xg) -> dataOut[0]
    dataOut.emplace_back(H.nobs(), 1);
    H.apply(xg, dataOut[0]);

    // H(Eg) -> dataOut[1]
    dataOut.emplace_back(H.nobs(), N);
    H.apply(Eg, dataOut[1]); 
}

void ESTKS::ensembleTransform(Ref<Array2d> E, vector<Array2d>& Egdata, 
                              const Ref<const Array> z, const CovarianceOperator& R, int k, 
                              Ref<Matrix> Xout) const
{
    double rho = 1.0 - (mInflation - 1.0);

    auto Hx = Egdata[0].matrix();
    auto HE = Egdata[1].matrix();

    int n = E.rows();
    int N = E.cols();
    int nobs = HE.rows();

    Matrix HL = HE * mT;  
    Matrix RinvHL(nobs, N-1);
    R.solve(HL, RinvHL); // R^ HL

    Matrix Ainv = HL.transpose() * RinvHL;   // N-1 x N-1
    Ainv.diagonal().array() += rho * (N-1);  // A^ = rho(N-1)I + (HL)' R^ HL

    
    // Compute the weight matrix `W` for updating ensemble perturbations

    //ENDAS_PERF_BEGIN(SymSqrt);
    Matrix C(Ainv.rows(), Ainv.cols());
    inverseSymmetricSqrt(Ainv, C, true);
    //ENDAS_PERF_END(SymSqrt);

    Matrix W = C * mT.transpose();
    W *= sqrt(N-1);

    // Compute the weight vector `w` for updating the ensemble mean

    auto dz = z - Hx.array();
    Matrix RinvDz(nobs, 1);
    R.solve(dz, RinvDz);

    ColVec w = HL.transpose() * RinvDz;     // N-1 x 1

    //ENDAS_PERF_BEGIN(AInvert);
    Eigen::LLT<Ref<Matrix>> AinvChol(Ainv);
    w = AinvChol.solve(w);
    //ENDAS_PERF_END(AInvert);

    W.colwise() += w;    

    Xout.array() = 1.0 / (double)N;
    Xout.noalias()+= mT * W;

    // Have the transform, apply it to the ensemble
    E.matrix() = E.matrix() * Xout;

    if (rho != 1.0)
    {
        Xout.array() = 1.0 / (double)N;
        Xout.noalias()+= (rho * mT) * W;
    }
}
