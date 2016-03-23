#include "Matrix.h"

Matrix::Matrix(int rows, int cols, Scalar val, Ordering ordering)
    :
      std::vector<Scalar>(rows*cols, val),
      ordering_(ordering),
      rows_(rows),
      cols_(cols)
{

}

Scalar& Matrix::operator ()(int row, int col)
{
    switch (ordering_)
    {
    case ROW_MAJOR: return std::vector<Scalar>::operator [](col + cols_*row);
    case COL_MAJOR: return std::vector<Scalar>::operator [](row + rows_*col);
    }

    throw std::exception();
}

const Scalar& Matrix::operator ()(int row, int col) const
{
    switch (ordering_)
    {
    case ROW_MAJOR: return std::vector<Scalar>::operator [](col + cols_*row);
    case COL_MAJOR: return std::vector<Scalar>::operator [](row + rows_*col);
    }

    throw std::exception();
}

void Matrix::resize(int rows, int cols)
{
    std::vector<Scalar>::resize(rows*cols);
    rows_ = rows;
    cols_ = cols;
}

Matrix& Matrix::operator+=(const Matrix& rhs)
{
    for(int i = 0, end = size(); i < end; ++i)
        operator [](i) += rhs[i];

    return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs)
{
    for(int i = 0, end = size(); i < end; ++i)
        operator [](i) -= rhs[i];

    return *this;
}

Matrix& Matrix::operator*=(Scalar rhs)
{
    for(int i = 0, end = size(); i < end; ++i)
        operator [](i) *= rhs;

    return *this;
}

Matrix& Matrix::operator/=(Scalar rhs)
{
    for(int i = 0, end = size(); i < end; ++i)
        operator [](i) /= rhs;

    return *this;
}

//- External functions

Matrix operator+(Matrix A, const Matrix& B)
{
    A += B;
    return A;
}

Matrix operator-(Matrix A, const Matrix& B)
{
    A -= B;
    return B;
}

Matrix operator*(Matrix A, Scalar a)
{
    A *= a;
    return A;
}

Matrix operator*(Scalar a, Matrix A)
{
    A *= a;
    return A;
}

