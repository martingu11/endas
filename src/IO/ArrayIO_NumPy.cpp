#include <Endas/IO/ArrayIO.hpp>
#include <Endas/Endas.hpp>
#include <cnpy.h>

using namespace std;
using namespace endas;


Array2d endas::loadArrayFromNpy(std::string path)
{
    cnpy::NpyArray npyA = cnpy::npy_load(path);

    ENDAS_REQUIRE(npyA.shape.size() <= 2, std::runtime_error, "Only 1 and 2 -dimensional arrays can be loaded");
    ENDAS_REQUIRE(npyA.word_size == sizeof(double), std::runtime_error, "Only double-precision real arrays can be loaded");

    int nrows = npyA.shape[0];
    int ncols = (npyA.shape.size() == 2)? npyA.shape[1] : 1;

    Array2d A(npyA.shape[0], npyA.shape[1]);

    const double* data = npyA.data<double>();

    // Already in fortran order -> memcpy
    if (npyA.fortran_order)
    {
        memcpy(A.data(), data, npyA.num_bytes());
    }
    // In C order -> copy transposed
    else
    {
        A = Eigen::Map<const Array2d>(data, nrows, ncols).transpose(); 
    }

    return A;
}

void endas::saveArrayAsNpy(const Ref<const Array2d> A, std::string path)
{
    vector<size_t> shape = { (size_t)A.rows(), (size_t)A.cols() };
    cnpy::npy_save(path, A.data(), shape, true);
}


