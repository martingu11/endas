#include <Endas/Utils/SecondOrderSampling.hpp>
#include <Endas/Core/Ensemble.hpp>
#include <Endas/Random/Random.hpp>

#include <Eigen/SVD>
#include <Eigen/QR>

#include <iostream>

using namespace std;
using namespace endas;


void endas::generateRandomOrthoMatrix(Ref<Array2d> out)
{
    ENDAS_ASSERT(out.rows() == out.cols());
    int n = out.rows();

    // Fill `out` with random numbers from N(0, 1)
    RandomNumberGenerator& rng = getRandomNumberGenerator();
    rng.standardNormal(out);

    // QR decomposition of `out` (in-place)
    Eigen::ColPivHouseholderQR<Ref<Matrix>> QR(out);

    Array x = QR.matrixR().matrix().diagonal();
    x /= x.cwiseAbs();

    out.matrix() = QR.householderQ();
    out.matrix() *= x.matrix().asDiagonal();
}




void SecondOrderExactSample::initFromStates(Ref<Array2d> states, bool subtractMean, double cutoff)
{
    if (subtractMean)
    {
        toAnomaly(states, states);
    }

    auto svd = Eigen::JacobiSVD<Matrix>(states.matrix(), Eigen::ComputeThinU);
    
    const auto& S = svd.singularValues();
    int numEOFs = svd.nonzeroSingularValues();

    // Using singular value cutoff. Look up the position of the singular value below the cutoff
    // threshold in the sorted `S` array. 
    if (cutoff > 0.0)
    {
        double maxS = S.maxCoeff(); 
        int a = 0;
        int b = numEOFs;
        while (a < b)
        {
            int mid = (a+b) / 2;
            auto s = S(mid) / maxS; // Normalize `S` to 0..1

            if (s > cutoff) a = mid + 1;
            else b = mid;
        }
        numEOFs = a;
    }

    auto n = states.rows();
    mS = svd.singularValues().block(0, 0, numEOFs, 1) * sqrt(1.0 / (double)(n - 1));
    mU = svd.matrixU().block(0, 0, n, numEOFs);
}

void SecondOrderExactSample::initFromEOFs(Array S, Array2d U)
{
    ENDAS_ASSERT(U.cols() == S.size());
    mS = move(S);
    mU = move(U);
}

const Array& SecondOrderExactSample::getS() const
{
    return mS;
}

const Array2d& SecondOrderExactSample::getU() const
{
    return mU;
}

index_t SecondOrderExactSample::numEOFs() const
{
    return mS.size();
}

void SecondOrderExactSample::samplePerturbations(Ref<Array2d> out)
{
    auto n = out.rows();
    auto N = out.cols();
    auto numEOFs = mS.size();

    ENDAS_CHECK_ARGUMENT(N <= numEOFs, "N cannot be larger than the number of EOFs");
    ENDAS_ASSERT(n == mU.rows());

    auto U = mU.block(0, 0, n, N-1).matrix();
    auto S = mS.block(0, 0, N-1, 1).matrix();

    // Generate random orthonormal matrix
    Matrix theta(N, N);
    generateRandomOrthoMatrix(theta);

    // Theta * Sigma where Sigma is a diagonal matrix with S on the diagonal
    Matrix thetaS = theta.block(0, 0, N, N-1) * S.asDiagonal();

    // Sampled perturbations as A' = (sqrt(N-1) * U * thetaS^T. This is a single DGEMM call 
    // if BLAS is enabled
    out.matrix().noalias() = sqrt(N-1) * U * thetaS.transpose();

    // Remove any mean from the perturbations
    toAnomaly(out, out);
}
