#include <Endas/Spatial/CoordinateSystem.hpp>
#include <Endas/Endas.hpp>

using namespace std;
using namespace endas;


constexpr double DEG2RAD = 0.0174532925; // Degree to radian conversion

constexpr Eigen::Index LAT = 0;
constexpr Eigen::Index LON = 1;


CoordinateSystem::~CoordinateSystem()
{ }


EuclideanCS::EuclideanCS(int ndim)
: mNdim(ndim)
{ }

int EuclideanCS::ndim() const 
{
    return mNdim;
}

bool EuclideanCS::isCartesian() const
{
    return true;
}

void EuclideanCS::distance(const Ref<const Array2d> A, const Ref<const Array2d> B,
                           Ref<Array> out) const
{
    ENDAS_ASSERT(A.cols() == B.cols());
    ENDAS_ASSERT(A.rows() == B.rows());
    ENDAS_ASSERT(A.cols() == out.size());

    out = (A - B).square().colwise().sum().sqrt();
}



LatLonCS::LatLonCS(double R)
: mR(R)
{ }

int LatLonCS::ndim() const 
{
    return 2;
}

bool LatLonCS::isCartesian() const
{
    return false;
}


inline double square(double a) { return a*a; }

inline double haversine(double Alat, double Alon, double Blat, double Blon, double R)
{
    Alat*= DEG2RAD;
    Blat*= DEG2RAD;
    Alon*= DEG2RAD;
    Blon*= DEG2RAD;

    double a = square(sin((Alat - Blat) / 2.0)) + cos(Alat) * cos(Blat) * square(sin((Alon - Blon) / 2.0));
    if (a > 1.0) a = 1.0;
    a = 2.0 * asin(sqrt(a));
    return a * R;
}


void LatLonCS::distance(const Ref<const Array2d> A, const Ref<const Array2d> B, Ref<Array> out) const
{
    ENDAS_ASSERT(A.cols() == B.cols());
    ENDAS_ASSERT(A.rows() == 2 && B.rows() == 2);
    ENDAS_ASSERT(A.cols() == out.size());

    for (int j = 0; j != A.cols(); j++)
    {
        auto Aj = A.col(j);
        auto Bj = B.col(j);
        out(j) = haversine(Aj(LAT), Aj(LON), Bj(LAT), Bj(LON), mR);
    }
}

