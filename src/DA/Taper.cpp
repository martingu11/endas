#include <Endas/DA/Taper.hpp>

using namespace std;
using namespace endas;


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
    ENDAS_NOT_IMPLEMENTED;
}



NoTaper::NoTaper(double L)
: TaperFn(L) 
{ }

void NoTaper::taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const
{
    ENDAS_NOT_IMPLEMENTED;
}


LinearTaper::LinearTaper(double L)
: TaperFn(L) 
{ }

void LinearTaper::taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const
{
    ENDAS_NOT_IMPLEMENTED;
}


SphericalTaper::SphericalTaper(double L)
: TaperFn(L) 
{ }

void SphericalTaper::taper(const Ref<const Array> x, const Ref<const Array> d, Ref<Array> out) const
{
    ENDAS_NOT_IMPLEMENTED;
}



