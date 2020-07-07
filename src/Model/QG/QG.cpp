#include <EndasModels/QG.hpp>
#include <Endas/Endas.hpp>
#include "../../Compatibility.hpp"

#include <iostream>


using namespace std;
using namespace endas;


static constexpr int QG_N = 129;
static constexpr int QG_M = 129;
static constexpr int QG_SIZE = QG_N*QG_M;

static bool globFmodInitialized = false;


struct QGModel::Data
{
    double rkb;
    double rkh;
    double rkh2;
    double F; 
    double r;

    double t;
    double internalDt;


    //Array2d Q;

    Data(int ensSize) 
    : rkb(0), 
      rkh(0),
      rkh2(2e-12),
      F(1600),
      r(1e-5),
      t(0)
    { 
        init();
    }

    void init();
};

QGModel::QGModel(int N, double internalStep)
: mData(make_unique<Data>(N))
{ 
    mData->internalDt = internalStep;
}



QGModel::~QGModel() { }

#if !ENDAS_HAS_FORTRAN

void QGModel::Data::init(int ensSize) 
{ 
    ENDAS_NOT_SUPPORTED("EnDAS must be compiled with Fortran support to enable the QG model.");
}

void QGModel::operator()(Ref<Array2d> x, int k, double dt, bool store) const
{ }   


#else

// Defined in qgstep.f90
extern "C" {

void qg_step_rk4(double t, double dt, double rkb, double rkh, double rkh2, 
                 double F, double r, double* PSI, double* Q);

void qg_params_init();

void qg_laplacian(const double* A, double dx, double dy, double* L);

void qg_calc_psi(const double* psiguess, const double* Q, double* psi, double F);

}

void QGModel::Data::init() 
{ 
    if (!globFmodInitialized)
    {
        qg_params_init();
        globFmodInitialized = true;
    }
}


int QGModel::sizex() const
{
    return QG_M;
}

int QGModel::sizey() const
{
    return QG_N;
}


void QGModel::operator()(Ref<Array2d> x, int k, double dt, bool store) const
{
    Array2d Q(QG_SIZE, x.cols());

    double dx = 1.0 / double(QG_N - 1);
    double dy = 1.0 / double(QG_M - 1);

    double t = mData->t;
    double tend = t + dt;

    #pragma omp parallel for
    for (int j = 0; j != x.cols(); j++)
    {
        qg_laplacian(x.col(j).data(), dx, dy, Q.col(j).data());
        Q.col(j) -= mData->F * x.col(j);

        double t = mData->t;
        double tend = t + dt;

        while (t < tend)
        {
            double thisdt = std::min(mData->internalDt, tend - t);

            qg_step_rk4(
                t, thisdt, mData->rkb, mData->rkh, mData->rkh2, mData->F, mData->r, 
                x.col(j).data(), Q.col(j).data());

            t+= thisdt;
        }

        Array2d xguess = x.col(j);
        qg_calc_psi(xguess.data(), Q.col(j).data(), x.col(j).data(), mData->F);
    }

    mData->t += dt;
}



#endif


