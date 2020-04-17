#include <Endas/Core/LinAlg.hpp>
#include <Endas/Endas.hpp>


using namespace std;
using namespace endas;

static const Matrix globEmptyMatrix;


const Matrix& endas::emptyMatrix()
{
    return globEmptyMatrix;
}



/*Vector endas::makeVector(std::initializer_list<real_t> values)
{
    Vector vec(values.size());

    real_t* ptr = vec.data();
    real_t* end = ptr + vec.size();
    for (auto&& x : values)
    {
        (*ptr++) = x;
    }
    return vec;
}*/


Array endas::makeArray(std::initializer_list<real_t> values)
{
    Array vec(values.size());

    real_t* ptr = vec.data();
    real_t* end = ptr + vec.size();
    for (auto&& x : values)
    {
        (*ptr++) = x;
    }
    return vec;
}


Matrix endas::makeMatrix(int rows, int cols, std::initializer_list<real_t> values)
{
    ENDAS_ASSERT(values.size() >= rows*cols);
    Matrix mat(rows, cols);

    // Not the most efficient perhaps but this is for small hand-written matrices anyway
    auto value = values.begin();
    for (int i = 0; i != mat.rows(); i++)
        for (int j = 0; j != mat.cols(); j++)
            mat(i, j) = *value++;
    return mat;
}



