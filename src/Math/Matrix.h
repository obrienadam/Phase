#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

#include "Types.h"

class Matrix : public std::vector<Scalar>
{
public:

    Matrix(Size m = 0, Size n = 0, const std::initializer_list<Scalar>& coeffs = {});

    void resize(Size m, Size n);

    void zero();

    void init(const Scalar *begin, const Scalar *end);

    Size m() const
    { return m_; }

    Size n() const
    { return n_; }

    bool isSquare() const
    { return m_ == n_; }

    void scaleRow(Size i, Scalar factor);

    Scalar &operator()(Size i, Size j);

    Scalar operator()(Size i, Size j) const;

    Matrix &operator=(const std::initializer_list<Scalar> &list);

    Matrix &operator+=(const Matrix &rhs);

    Matrix &operator-=(const Matrix &rhs);

    Matrix &operator*=(Scalar rhs);

    Matrix &operator/=(Scalar rhs);

    Matrix &solve(Matrix &b);

    Matrix &transpose();

    Matrix &invert();

    Scalar norm(char type = 'I') const;

    Scalar cond(char type = 'I') const;

    Matrix subMatrix(size_t startRow, size_t startCol, size_t endRow, size_t endCol) const;

    Scalar *data()
    { return std::vector<Scalar>::data(); }

    const Scalar *data() const
    { return std::vector<Scalar>::data(); }

    std::vector<Scalar> containerCopy() const
    { return std::vector<Scalar>(*this); }

private:

    Size m_, n_;
    mutable std::vector<int> ipiv_;

};

//- External functions
Matrix eye(size_t nRows, size_t nCols);

Matrix random(size_t nRows, size_t nCols, Scalar min, Scalar max);

Matrix transpose(Matrix mat);

Matrix inverse(Matrix mat);

Matrix solve(Matrix A, Matrix b);

Matrix sum(const Matrix &A);

Matrix operator+(Matrix lhs, const Matrix &rhs);

Matrix operator-(Matrix lhs, const Matrix &rhs);

Matrix operator*(Matrix lhs, Scalar rhs);

Matrix operator*(Scalar lhs, Matrix rhs);

Matrix operator*(const Matrix &lhs, const Matrix &rhs);

Matrix operator/(Matrix lhs, Scalar rhs);

std::ostream &operator<<(std::ostream &os, const Matrix &mat);

#endif
