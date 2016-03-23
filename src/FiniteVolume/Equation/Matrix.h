#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

#include "Types.h"

class Matrix : protected std::vector<Scalar>
{
public:

    enum Ordering{ROW_MAJOR, COL_MAJOR};

    Matrix(Ordering ordering = ROW_MAJOR) : ordering_(ordering) {}
    Matrix(int rows, int cols, Scalar val = 0., Ordering ordering = ROW_MAJOR);

    int rows() const { return rows_; }
    int cols() const { return cols_; }

    Scalar& operator()(int row, int col);
    const Scalar& operator()(int row, int col) const;

    void resize(int rows, int cols);

    Matrix& operator+=(const Matrix& rhs);
    Matrix& operator-=(const Matrix& rhs);
    Matrix& operator*=(Scalar rhs);
    Matrix& operator/=(Scalar rhs);

private:

    Ordering ordering_;
    int rows_, cols_;

};

Matrix operator+(Matrix A, const Matrix& B);
Matrix operator-(Matrix A, const Matrix& B);
Matrix operator*(Matrix A, Scalar a);
Matrix operator*(Scalar a, Matrix A);

#endif
