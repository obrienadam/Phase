#include <ostream>

#ifdef __INTEL_COMPILER
#include <mkl.h>
#else
extern "C"
{
#define lapack_complex_float float _Complex
#define lapack_complex_double double _Complex
#include <lapacke.h>
#include <cblas.h>
}
#endif

#include "Matrix.h"
#include "Exception.h"

Matrix::Matrix(size_t nRows, size_t nCols)
{
    resize(nRows, nCols);
}

void Matrix::resize(size_t nRows, size_t nCols)
{
    nRows_ = nRows;
    nCols_ = nCols;
    isSquare_ = nRows_ == nCols_;
    std::vector<Scalar>::resize(nRows*nCols);
    ipiv_.resize(nRows_);
    zero();
}

void Matrix::zero()
{
    std::fill(begin(), end(), 0.);
}

void Matrix::init(const Scalar *begin, const Scalar *end)
{
    std::vector<Scalar>::assign(begin, end);
}

Scalar& Matrix::operator()(size_t m, size_t n)
{
    return std::vector<Scalar>::operator [](m*nCols_ + n);
}

const Scalar& Matrix::operator()(size_t m, size_t n) const
{
    return std::vector<Scalar>::operator [](m*nCols_ + n);
}

Matrix& Matrix::operator=(const std::initializer_list<Scalar>& list)
{
    if(list.size() != nRows_*nCols_)
        throw Exception("Matrix", "operator=", "initializer list size does not match matrix dimensions.");

    std::vector<Scalar>::operator =(list);
    return *this;
}

//- Operators
Matrix& Matrix::operator+=(const Matrix& rhs)
{
    for(size_t m = 0; m < nRows_; ++m)
        for(size_t n = 0; n < nCols_; ++n)
            operator ()(m, n) += rhs(m, n);

    return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs)
{
    for(size_t m = 0; m < nRows_; ++m)
        for(size_t n = 0; n < nCols_; ++n)
            operator ()(m, n) -= rhs(m, n);

    return *this;
}

Matrix& Matrix::operator*=(Scalar rhs)
{
    for(size_t m = 0; m < nRows_; ++m)
        for(size_t n = 0; n < nCols_; ++n)
            operator ()(m, n) *= rhs;

    return *this;
}

Matrix& Matrix::operator/=(Scalar rhs)
{
    for(size_t m = 0; m < nRows_; ++m)
        for(size_t n = 0; n < nCols_; ++n)
            operator ()(m, n) /= rhs;

    return *this;
}

Matrix& Matrix::solve(Matrix& b)
{
    if(isSquare_)
        LAPACKE_dgesv(LAPACK_ROW_MAJOR, nRows_, b.nCols_, data(), nCols_, ipiv_.data(), b.data(), b.nCols_);
    else
    {
        LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', nRows_, nCols_, b.nCols_, data(), nCols_, b.data(), b.nCols_);
        b.nRows_ = nCols_; // A bit hackish, but resizes b apropriately. Only works so long as the data format is row major
    }

    return b;
}

Matrix& Matrix::transpose()
{
    if(isSquare_) // Square matrices
    {
        auto &self = *this;

        for(size_t m = 0; m < nRows_; ++m)
            for(size_t n = m + 1; n < nCols_; ++n)
                std::swap(self(m, n), self(n, m));
    }
    else // General non-square matrices
    {
        auto last = end();
        auto first = begin();
        int m = nCols_;

        const int mn1 = (last - first - 1);
        const int n   = (last - first) / m;
        std::vector<bool> visited(last - first);
        auto cycle = first;
        while (++cycle != last)
        {
            if (visited[cycle - first])
                continue;
            int a = cycle - first;
            do  {
                a = a == mn1 ? mn1 : (n * a) % mn1;
                std::swap(*(first + a), *cycle);
                visited[a] = true;
            } while ((first + a) != cycle);
        }

        std::swap(nRows_, nCols_);
    }

    return *this;
}

Matrix& Matrix::invert()
{
    if(nRows_ != nCols_)
        operator=(Matrix(*this).transpose()*(*this));

    lapack_int info1 = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nRows_, nCols_, data(), nCols_, ipiv_.data());
    lapack_int info2 = LAPACKE_dgetri(LAPACK_ROW_MAJOR, nRows_, data(), nCols_, ipiv_.data());

    if(info1 != 0 || info2 != 0)
        throw Exception("Matrix", "invert", "inversion failed, matrix is singular to working precision.");

    return *this;
}

Scalar Matrix::norm(char type) const
{
    return LAPACKE_dlange(LAPACK_ROW_MAJOR, type, nRows_, nCols_, data(), nCols_);
}

Scalar Matrix::cond(char type) const
{
    Scalar rcond;
    LAPACKE_dgecon(LAPACK_ROW_MAJOR, type, nRows_, data(), nCols_, norm(type), &rcond);
    return rcond;
}

Matrix Matrix::subMatrix(size_t startRow, size_t startCol, size_t endRow, size_t endCol) const
{
    Matrix subMat(endRow - startRow, endCol - startCol);
    auto &self = *this;

    for(size_t m = startRow; m < endRow; ++m)
        for(size_t n = startCol; n < endCol; ++n)
            subMat(m - startRow, n - startCol) = self(m, n);

    return subMat;
}

//- External functions
Matrix eye(size_t nRows, size_t nCols)
{
    Matrix mat(nRows, nCols);
    for(size_t i = 0; i < nRows && i < nCols; ++i)
        mat(i, i) = 1.;

    return mat;
}

Matrix random(size_t nRows, size_t nCols, Scalar min, Scalar max)
{
    Matrix mat(nRows, nCols);

    for(size_t m = 0; m < nRows; ++m)
        for(size_t n = 0; n < nCols; ++n)
            mat(m, n) = (double)rand()/RAND_MAX*(max - min) + min;

    return mat;
}

Matrix transpose(Matrix mat)
{
    return mat.transpose();
}

Matrix inverse(Matrix mat)
{
    return mat.invert();
}

Matrix solve(Matrix A, Matrix b)
{
    return A.solve(b);
}

Matrix sum(const Matrix& A)
{
    Matrix result(1, A.nRows() == 1 ? 1 : A.nCols());
    result.zero();

    if(result.nCols() == 1)
    {
        for(int i = 0; i < A.nCols(); ++i)
            result(0, 0) += A(0, i);
    }
    else
    {
        for(int i = 0; i < A.nRows(); ++i)
            for(int j = 0; j < A.nCols(); ++j)
                result(0, j) += A(i, j);
    }

    return result;
}

Matrix operator+(Matrix lhs, const Matrix& rhs)
{
    return lhs += rhs;
}

Matrix operator-(Matrix lhs, const Matrix& rhs)
{
    return lhs -= rhs;
}

Matrix operator*(Matrix lhs, Scalar rhs)
{
    return lhs *= rhs;
}

Matrix operator*(Scalar lhs, Matrix rhs)
{
    return rhs *= lhs;
}

Matrix operator*(const Matrix& A, const Matrix& B)
{
    Matrix C(A.nRows(), B.nCols());

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A.nRows(), B.nCols(), A.nCols(), 1., A.data(), A.nCols(), B.data(), B.nCols(), 1., C.data(), C.nCols());

    // Works for sure
    //    const int nI = A.nRows(), nJ = B.nCols(), nK = A.nCols();

    //    for(int i = 0; i < nI; ++i)
    //        for(int j = 0; j < nJ; ++j)
    //            for(int k = 0; k < nK; ++k)
    //                C(i, j) += A(i, k)*B(k, j);

    return C;
}

Matrix operator/(Matrix lhs, Scalar rhs)
{
    return lhs /= rhs;
}

std::ostream& operator<<(std::ostream& os, const Matrix& mat)
{
    for(size_t m = 0; m < mat.nRows(); ++m)
    {
        for(size_t n = 0; n < mat.nCols(); ++n)
        {
            os << mat(m, n) << ' ';
        }
        os << '\n';
    }

    return os;
}
