#ifndef PHASE_MATRIX_H
#define PHASE_MATRIX_H

#include <vector>

#include "Types/Types.h"

class Matrix
{
public:

    Matrix(Size m = 0, Size n = 0, const std::initializer_list<Scalar> &coeffs = {});

    void resize(Size m, Size n);

    void zero();

    void setIdentity(Size m);

    void init(const Scalar *begin, const Scalar *end);

    Size m() const
    { return m_; }

    Size n() const
    { return n_; }

    bool isSquare() const
    { return m_ == n_; }

    void setRow(int i, const std::initializer_list<Scalar> &coeffs);

    void addToRow(int i, const std::initializer_list<Scalar> &coeffs);

    void scaleRow(Size i, Scalar factor);

    Scalar &operator()(Size i, Size j);

    Scalar operator()(Size i, Size j) const;

    Matrix &operator=(const std::initializer_list<Scalar> &list);

    Matrix &operator+=(const Matrix &rhs);

    Matrix &operator-=(const Matrix &rhs);

    Matrix &operator*=(Scalar rhs);

    Matrix &operator/=(Scalar rhs);

    Matrix &solve(Matrix &b);

    Matrix solve(const Matrix &b);

    Matrix &transpose();

    Matrix &invert();

    Matrix &pinvert();

    Scalar norm(char type = 'I') const;

    Scalar cond(char type = 'I') const;

    Matrix subMatrix(size_t startRow, size_t startCol, size_t endRow, size_t endCol) const;

    Scalar *data()
    { return vals_.data(); }

    const Scalar *data() const
    { return vals_.data(); }

    std::vector<Scalar>::const_iterator begin() const
    { return vals_.begin(); }

    std::vector<Scalar>::const_iterator end() const
    { return vals_.end(); }

private:

    static Matrix _tmp;

    Size m_, n_;

    std::vector<Scalar> vals_;

    mutable std::vector<int> ipiv_;

};

//- External functions
Matrix eye(int m);

Matrix random(size_t nRows, size_t nCols, Scalar min, Scalar max);

Matrix transpose(const Matrix &mat);

Matrix inverse(Matrix mat);

Matrix pseudoInverse(Matrix mat);

Matrix solve(Matrix A, Matrix b);

Matrix sum(const Matrix &A);

Matrix operator+(Matrix lhs, const Matrix &rhs);

Matrix operator-(Matrix lhs, const Matrix &rhs);

Matrix operator*(Matrix lhs, Scalar rhs);

Matrix operator*(Scalar lhs, Matrix rhs);

Matrix operator*(const Matrix &lhs, const Matrix &rhs);

Matrix operator/(Matrix lhs, Scalar rhs);

Matrix multiply(const Matrix &A, const Matrix &B, bool transA = false, bool transB = false);

std::ostream &operator<<(std::ostream &os, const Matrix &mat);

#endif
