#include <Endas/DA/Taper.hpp>

using namespace std;
using namespace endas;

inline double pow2(double x) { return x*x; }
inline double pow3(double x) { return x*x*x; }
inline double pow4(double x) { return x*x*x*x; }
inline double pow5(double x) { return x*x*x*x*x; }



TaperFn::TaperFn(double L)
: mL(L)
{ }


TaperFn::~TaperFn()
{ }

double TaperFn::supportRange() const { return mL; }


GaspariCohnTaper::GaspariCohnTaper(double L) 
: TaperFn(L) 
{ }

void GaspariCohnTaper::taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const
{
    double L = this->mL;

    static constexpr double FRAC_5_3 = 5/3.0;
    static constexpr double FRAC_5_8 = 5/8.0;
    static constexpr double FRAC_1_2 = 0.5;
    static constexpr double FRAC_1_4 = 0.25;
    static constexpr double FRAC_1_12 = 1/12.0;

    out = x.binaryExpr(d, [=](real_t x, real_t d) -> real_t
    {
        double r = d / L;
        if (r < 1) 
            return x * (1.0         - FRAC_5_3*pow2(r) + FRAC_5_8*pow3(r) + FRAC_1_2*pow4(r) - FRAC_1_4*pow5(r));
        else if (r < 2)
            return x * (4.0 - 5.0*r + FRAC_5_3*pow2(r) + FRAC_5_8*pow3(r) - FRAC_1_2*pow4(r) + FRAC_1_12*pow5(r) - 2.0/(3.0*r));
        else
            return 0;
    });
}


NoTaper::NoTaper(double L)
: TaperFn(L) 
{ }

void NoTaper::taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const
{
    if (out.data() != x.data())
    {
        out = x;
    }
}


LinearTaper::LinearTaper(double L)
: TaperFn(L) 
{ }

void LinearTaper::taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const
{
    double L = this->mL;
    out = x.binaryExpr(d, [=](real_t x, real_t d) -> real_t
    {
        double r = d / L;
        return (r < 1)? x * (1.0 - r) : 0;
    });
}


SphericalTaper::SphericalTaper(double L)
: TaperFn(L) 
{ }

void SphericalTaper::taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const
{
    double L = this->mL;
    out = x.binaryExpr(d, [=](real_t x, real_t d) -> real_t
    {
        double r = d / L;
        return (r < 1)? x * (1.0 - (1.5*r - 0.5*pow3(r))) : 0;
    });
}



