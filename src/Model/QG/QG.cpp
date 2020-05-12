#include <EndasModels/QG.hpp>
#include <Endas/Endas.hpp>
#include "../../Compatibility.hpp"

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

    Array2d Q;

    Data(int ensSize) 
    : rkb(0), 
      rkh(0),
      rkh2(2e-12),
      F(1600),
      r(1e-5),
      t(0)
    { 
        
    }

    void init(const Ref<const Array2d> x0);
};

QGModel::QGModel(int N)
: mData(make_unique<Data>(N))
{ 
}



QGModel::~QGModel() { }

#if !ENDAS_HAS_FORTRAN

void QGModel::Data::init(int ensSize) { }

void QGModel::operator()(Ref<Array2d> x, int k, double dt, bool store) const
{
    ENDAS_NOT_SUPPORTED("The QG model has not been compiled into EnDAS.");
}

#else

// Defined in qgstep.f90
extern "C" {

void qg_step_rk4(double t, double dt, double rkb, double rkh, double rkh2, 
                 double F, double r, double* PSI, double* Q);

void qg_params_init();

void qg_laplacian(const double* A, double dx, double dy, double* L);

void qg_calc_psi(const double* psiguess, const double* Q, double* psi, double F);

}


void QGModel::Data::init(const Ref<const Array2d> x0) 
{ 
    if (!globFmodInitialized)
    {
        qg_params_init();
        globFmodInitialized = true;
    }

    Q = Array2d(QG_SIZE, x0.cols());


    double dx = 1.0 / double(QG_N - 1);
    double dy = 1.0 / double(QG_M - 1);

    for (int j = 0; j != x0.cols(); j++)
    {
        qg_laplacian(x0.col(j).data(), dx, dy, Q.col(j).data());
        Q.col(j) -= F * x0.col(j);
    }


}


void QGModel::init(const Ref<const Array2d> x0)
{
    mData->init(x0);
}



void QGModel::operator()(Ref<Array2d> x, int k, double dt, bool store) const
{
    for (int j = 0; j != x.cols(); j++)
    {
        qg_step_rk4(
            mData->t, dt, mData->rkb, mData->rkh, mData->rkh2, mData->F, mData->r, 
            x.col(j).data(), mData->Q.col(j).data());
    };

    mData->t += dt;
}


void QGModel::calc_psi(const Ref<const Array2d> E, Ref<Array2d> Eout)
{
    for (int j = 0; j != E.cols(); j++)
    {
        qg_calc_psi(E.col(j).data(), mData->Q.col(j).data(), Eout.col(j).data(), mData->F);
    }
}



#endif


