#include <Endas/Core/LinAlg.hpp>
#include <Endas/Endas.hpp>

#include <Eigen/SVD>

using namespace std;
using namespace endas;



const Array& endas::emptyArray()
{
    static const Array emptyArray;
    return emptyArray;
}

const Array2d& endas::emptyArray2d()
{
    static const Array2d emptyArray2d;
    return emptyArray2d;
}


const Matrix& endas::emptyMatrix()
{
    static const Matrix emptyMatrix;
    return emptyMatrix;
}


Array endas::makeSequence(real_t start, real_t end, real_t step)
{
    index_t size = (index_t)ceil((end - start) / step);

    Array A(size);
    index_t i = 0;
    while (start < end)
    {
        A(i++) = start;
        start+= step;
    }
    return move(A);
}




Array endas::makeArray(std::initializer_list<real_t> values)
{
    Array vec(values.size());

    real_t* ptr = vec.data();
    real_t* end = ptr + vec.size();
    for (auto&& x : values)
    {
        (*ptr++) = x;
    }
    return move(vec);
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
    
    return move(mat);
}

ArrayShape endas::makeShape(index_t a)
{
    ArrayShape shape(1);
    shape << a;
    return shape;
}

ArrayShape endas::makeShape(index_t a, index_t b)
{
    ArrayShape shape(2);
    shape << a, b;
    return shape;
}

ArrayShape endas::makeShape(index_t a, index_t b, index_t c)
{
    ArrayShape shape(3);
    shape << a, b, c;
    return shape;
}

ArrayShape endas::makeShape(index_t a, index_t b, index_t c, index_t d)
{
    ArrayShape shape(4);
    shape << a, b, c, d;
    return shape;
}

ArrayShape2d endas::makeShape2d(index_t a, index_t b)
{
    ArrayShape2d shape;
    shape << a, b;
    return shape;
}

ArrayShape3d endas::makeShape3d(index_t a, index_t b, index_t c)
{
    ArrayShape3d shape;
    shape << a, b, c;
    return shape;
}




void endas::select(const Ref<const Array> A, const IndexArray& indices, Ref<Array> out)
{
    index_t ii = 0;
    for (index_t i : indices) out(ii++) = A(i);
}

void endas::selectRows(const Ref<const Array2d> A, const IndexArray& rows, Ref<Array2d> out)
{
    index_t rr = 0;
    for (index_t r : rows) out.row(rr++) = A.row(r);
}

void endas::distributeRows(const Ref<const Array2d> A, const IndexArray& rows, Ref<Array2d> out)
{
    index_t rr = 0;
    for (index_t r : rows) out.row(r) = A.row(rr++);
}


void endas::selectCols(const Ref<const Array2d> A, const IndexArray& cols, Ref<Array2d> out)
{
    index_t cc = 0;
    for (index_t c : cols) out.col(cc++) = A.col(c);

}

void endas::distributeCols(const Ref<const Array2d> A, const IndexArray& cols, Ref<Array2d> out)
{
    index_t cc = 0;
    for (index_t c : cols) out.col(c) = A.col(cc++);
}


void endas::selectRowsCols(const Ref<const Array2d> A, const IndexArray& rows, 
                           const IndexArray& cols, Ref<Array2d> out)
{
    index_t ii = 0;
    for (index_t i : cols)
    {
        index_t jj = 0;
        for (index_t j : rows)
        {
            out(ii, jj++) = A(i, j);
        }
        ++ii;
    }
}



void endas::inverseSymmetricSqrt(const Ref<const Matrix> A, Ref<Matrix> out, bool noalias)
{
    auto Asvd = A.jacobiSvd(Eigen::ComputeFullU);

    const Matrix& U = Asvd.matrixU();
    Array s = Asvd.singularValues().array().pow(-0.5);

    Matrix SU = s.matrix().asDiagonal() * U.transpose();

    if (noalias)
    {
        out.noalias() = U * SU;
    }
    else
    {
        out = U * SU;        
    }
}


