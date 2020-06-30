
#include <Endas/Random/Random.hpp>
#include <Endas/Random/GaussianRandomField.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"

#include <simple_fft/fft.hpp>

#include <iostream>

using namespace std;
using namespace endas;

static constexpr int MAX_SIZE = 8192*8192;


struct GaussianRandomField::Data
{
    int nx, ny;
    int s;

    Array2d eigenValues;
    double eigMinCoeff;

    void init(int nx, int ny, const IsotropicCovarianceFn& covFn, double eigMin);

};


void GaussianRandomField::Data::init(int nx, int ny, const IsotropicCovarianceFn& covFn, double eigMin)
{
    // Bump nx and ny up to nearest power of two

    int nybase = (int)ceil(log2(ny));
    int nxbase = (int)ceil(log2(nx));

    ny = pow(2, nybase);
    nx = pow(2, nxbase);

    cout << ny << "  " << nx << endl;
    

    if (nx*ny > MAX_SIZE)
        throw runtime_error("Field dimensions too high for GaussianRandomField");

    for (int s = 1; s < 10; s++)
    {
        if (nx*ny > MAX_SIZE) break;

        Array tx = makeSequence(0, nx);
        Array ty = makeSequence(0, ny);

        Array2d cov = Array2d::Zero(ny, nx);

        // Sample the covariance function at grid points, relative to 0,0 
        /// @todo Here we could avoid redoing what we already did in earlier loops

        double ty0 = ty(0);
        for (int j = 0; j != nx; j++)
        {
            double dtx = tx(j) - tx(0);
            double dtxSq = dtx * dtx;
            auto h = ty.unaryExpr([=](real_t _ty) 
            { 
                double dty = _ty - ty0;
                return sqrt(dtxSq + (dty*dty));
            });

            covFn.values(h, cov.col(j));
        }

        ComplexArray2d C_BCCB(ny*2, nx*2);
        
        C_BCCB.imag() = 0;
        C_BCCB.block(0, 0, ny, nx).real() = cov;
        C_BCCB.block(0, nx, ny, nx).real() = cov.rowwise().reverse(); // Horizontal flip
        C_BCCB.block(ny, 0, ny, 2*nx) = C_BCCB.block(0, 0, ny, 2*nx).colwise().reverse();

        
        // In-place 2d FFT of C_BCCB
        const char* error;
        if (!simple_fft::FFT(C_BCCB, ny*2, nx*2, error))
        {
            throw runtime_error(error);
        }

        auto eigVals = C_BCCB.real() / (2*ny-1) / (2*nx-1);

        // Check if eigenvalues are positive. If yes (or close to zero), we're fine. Otherwise
        // increase nx and ny and try again...
        this->eigMinCoeff = eigVals.minCoeff();

        this->nx = nx;
        this->ny = ny;
        this->s = s;

        if (this->eigMinCoeff >= eigMin)
        {
            // Truncate any remaining negative eigenvalues
            this->eigenValues = eigVals.unaryExpr([](real_t x)
            {
                return sqrt(std::max(0.0, x));
            });
            // We're done
            return;
        }

        nx*= 2;
        ny*= 2;
    }    

    throw runtime_error("Cannot instantiate GaussianRandomField: circulant embedding failed");

}



GaussianRandomField::GaussianRandomField(int nx, int ny, const IsotropicCovarianceFn& covFn, double eigMin)
: mData(make_unique<Data>())
{ 
    mData->init(nx, ny, covFn, eigMin);
}

GaussianRandomField::~GaussianRandomField()
{ }

void GaussianRandomField::sample(Ref<Array2d> out) const
{
    Array2d dummy;
    return this->sample(out, dummy);
}

void GaussianRandomField::sample(Ref<Array2d> out, Ref<Array2d> out2) const
{
    int nx = mData->nx;
    int ny = mData->ny;

    RandomNumberGenerator& rng = getRandomNumberGenerator();
    
    ComplexArray2d A(ny*2, nx*2);

    //ENDAS_PERF_BEGIN(GRFrand);
    rng.standardNormal(A.real());
    rng.standardNormal(A.imag());
    //ENDAS_PERF_END(GRFrand);

    // Generate the random Gausian field, we get two independent realizations in .real and .imag
    A*= mData->eigenValues;

    // In-place 2d FFT of C_BCCB
    //ENDAS_PERF_BEGIN(GRFFFT);
    const char* error;
    if (!simple_fft::FFT(A, ny*2, nx*2, error))
    {
        throw runtime_error(error);
    }
    //ENDAS_PERF_END(GRFFFT);

    auto Aout = A.block(0, 0, out.rows(), out.cols());
    out = Aout.real();
    if (out2.size() > 0) out2 = Aout.imag();
}