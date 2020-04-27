#include <Endas/Spatial/Variogram.hpp>
#include <Endas/Spatial/CoordinateSystem.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;


CovarianceFn::~CovarianceFn()
{ }

IsotropicCovarianceFn::~IsotropicCovarianceFn()
{ }


void IsotropicCovarianceFn::values(const Ref<const Array2d> A, const Ref<const Array2d> B, 
                                   Ref<Array> out) const
{
    EuclideanCS cs(A.rows());
    cs.distance(A, B, out);
    this->values(out, out);
}



//---------------------------------------------------------------------------------------
// ExponentialFamilyCovFn
//---------------------------------------------------------------------------------------

ExponentialFamilyCovFn::ExponentialFamilyCovFn(double alpha, double L, double sigma)
: mAlpha(alpha), mL(L), mSigma(sigma)
{ }


void ExponentialFamilyCovFn::values(const Ref<const Array> h, Ref<Array> out) const
{
    double alpha = mAlpha;
    double L = mL;
    double sigma = mSigma;

    if (alpha == 1.0)
    {
        out = h.unaryExpr([=](real_t h) 
        { 
            return sigma * exp(-h/L); 
        });
    }
    else if (alpha == 2.0)
    {
        out = h.unaryExpr([=](real_t h) 
        { 
            double hDivL = h/L;
            return sigma * exp(-hDivL*hDivL); 
        });
    }
    else
    {
        out = h.unaryExpr([=](real_t h) 
        { 
            return sigma * exp(-pow(h/L, alpha)); 
        });
    }
}


//---------------------------------------------------------------------------------------
// SphericalCovFn
//---------------------------------------------------------------------------------------

SphericalCovFn::SphericalCovFn(double L, double sigma)
: mL(L), mSigma(sigma)
{ }


void SphericalCovFn::values(const Ref<const Array> h, Ref<Array> out) const
{
    double L = mL;
    double sigma = mSigma;

    out = h.unaryExpr([=](real_t h) 
    { 
        if (h < L)
        {
            double hDivL = h/L;
            return sigma * (1.0 - (1.5*hDivL - 0.5*pow(hDivL, 3)));
        }
        else
            return 0.0;
    });
}