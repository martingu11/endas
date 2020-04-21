#include <Endas/Model/Lorenz95.hpp>
#include <Endas/Endas.hpp>
#include "../Compatibility.hpp"

#include <unordered_map>
#include <iostream>

using namespace std;
using namespace endas;


// Data needed for linearization at x.
struct L95Trajectory
{
    double dt;
    Array2d x;
    Array2d k1, k2, k3, k4;

    L95Trajectory(int n, int N, double _dt)
    : x(n, N), k1(n, N), k2(n, N), k3(n, N), k4(n, N), dt(_dt)
    { }

    L95Trajectory(const L95Trajectory&) = delete;
    L95Trajectory& operator=(const L95Trajectory&) = delete;

    L95Trajectory(L95Trajectory&&) = default;
    L95Trajectory& operator=(L95Trajectory&&) = default;
};


struct Lorenz95Model::Data
{
    int n;
    int F;

    unordered_map<int, L95Trajectory> trajectories;

    int index(int i) { return (i < 0)? n+i : ((i >= n)? i-n : i); }

    template <class OutT, class XT>
    void l95(Eigen::ArrayBase<OutT>& out, const Eigen::ArrayBase<XT>& x, double dt)
    {
        for (int i = 0; i != n; i++)
        {
            int im2 = index(i-2);
            int im1 = index(i-1);
            int ip1 = index(i+1);

            out(i) = -x(im2)*x(im1) + x(im1)*x(ip1) - x(i) + F;
            out(i)*= dt;
        }
    }

    template <class OutT, class XT, class DXT>
    void l95tl(Eigen::ArrayBase<OutT>& out, const Eigen::ArrayBase<XT>& x, 
               const Eigen::ArrayBase<DXT>& dx, double dt)
    {
        for (int i = 0; i != n; i++)
        {
            int im2 = index(i-2);
            int im1 = index(i-1);
            int ip1 = index(i+1);
            int ip2 = index(i+2);

            out(i) = -x(im1)*dx(im2) + (x(ip1)-x(im2))*dx(im1) - dx(i) + x(im1)*dx(ip1);
            out(i)*= dt;
        }
    }

    template <class XT, class DXT>
    Array l95ad(const Eigen::ArrayBase<XT>& x, const Eigen::ArrayBase<DXT>& dx)
    {
        Array out(x.size());
        for (int i = 0; i != n; i++)
        {
            int im2 = index(i-2);
            int im1 = index(i-1);
            int ip1 = index(i+1);
            int ip2 = index(i+2);

            out(i) = x(im2)*dx(im1) + (x(ip2)-x(im1))*dx(ip1) - dx(i) - x(ip1)*dx(ip2);
        }
        return out;
    }



};


Lorenz95Model::Lorenz95Model(int n, int F)
: mData(make_unique<Data>())
{ 
    mData->n = n;
    mData->F = F;
}

Lorenz95Model::~Lorenz95Model() { }


void Lorenz95Model::apply(Ref<Array2d> x, int k, double dt, bool store) const
{
    int n = x.rows();
    int N = x.cols();
    ENDAS_ASSERT(n == mData->n);

    // State propagation using the 4th-order Runge-Kutta method. Store what we need for the 
    // tl() and adj() computations
    L95Trajectory trj(n, N, dt);
   
    for (int i = 0; i != N; i++)
    {
        auto xi = x.col(i);
        auto k1i = trj.k1.col(i);
        auto k2i = trj.k2.col(i);
        auto k3i = trj.k3.col(i);
        auto k4i = trj.k4.col(i);

        mData->l95(k1i, xi, dt);
        mData->l95(k2i, xi + k1i / 2.0, dt);
        mData->l95(k3i, xi + k2i / 2.0, dt);
        mData->l95(k4i, xi + k3i, dt);

        xi += (k1i + k2i*2.0 + k3i*2.0 + k4i) / 6.0;
        trj.x.col(i) = xi;
    }

    if (store)
    {
        /// @todo Not perfect although this is a toy model anyway.
        mData->trajectories.erase(k);
        mData->trajectories.emplace(k, move(trj));
    }
}

void Lorenz95Model::tl(Ref<Array2d> x, int k) const
{
    int n = x.rows();
    int N = x.cols();
    ENDAS_ASSERT(n == mData->n);

    auto it = mData->trajectories.find(k);
    ENDAS_ASSERT(it != mData->trajectories.end());

    const L95Trajectory& trj = it->second;

    Array dK1(n), dK2(n), dK3(n), dK4(n);

    for (int i = 0; i != N; i++)
    {
        int ti = i % trj.x.cols();
        auto xi = x.col(i);
        auto txi = trj.x.col(ti);
        auto k1i = trj.k1.col(ti);
        auto k2i = trj.k2.col(ti);
        auto k3i = trj.k3.col(ti);
        auto k4i = trj.k4.col(ti);

        mData->l95tl(dK1, txi, xi, trj.dt);
        mData->l95tl(dK2, txi + k1i / 2.0, xi + dK1 / 2.0, trj.dt);
        mData->l95tl(dK3, txi + k2i / 2.0, xi + dK2 / 2.0, trj.dt);
        mData->l95tl(dK4, txi + k3i, xi + dK3, trj.dt);

        xi += (dK1 + dK2*2.0 + dK3*2.0 + dK4) / 6.0;
    }
}

void Lorenz95Model::adj(Ref<Array2d> x, int k) const
{
    int n = x.rows();
    int N = x.cols();
    ENDAS_ASSERT(n == mData->n);

    auto it = mData->trajectories.find(k);
    ENDAS_ASSERT(it != mData->trajectories.end());

    const L95Trajectory& trj = it->second;

    auto adK1 = [&](int i, const Ref<const Array> dx) -> Array
    {
        return mData->l95ad(trj.x.col(i), dx) * trj.dt;
    };

    auto adK2 = [&](int i, const Ref<const Array> dx) -> Array
    {
        Array aux = mData->l95ad(trj.x.col(i) + trj.k1.col(i) / 2.0, dx);
        return (aux + adK1(i, aux) / 2.0) * trj.dt;
    };

    auto adK3 = [&](int i, const Ref<const Array> dx) -> Array
    {
        Array aux = mData->l95ad(trj.x.col(i) + trj.k2.col(i) / 2.0, dx);
        return (aux + adK2(i, aux) / 2.0) * trj.dt;
    };

    auto adK4 = [&](int i, const Ref<const Array> dx) -> Array
    {
        Array aux = mData->l95ad(trj.x.col(i) + trj.k3.col(i), dx);
        return (aux + adK3(i, aux)) * trj.dt;
    };

    for (int i = 0; i != N; i++)
    {
        int ti = i % trj.x.cols();
        auto xi = x.col(i);
        xi+= (adK1(ti, xi) + adK2(ti, xi)*2.0 + adK3(ti, xi)*2.0 + adK4(ti, xi)) / 6.0;
    }

}


void Lorenz95Model::stepFinished(int k) const 
{ 
    mData->trajectories.erase(k);
}



