#ifndef STATIC_MATRIX_H
#define STATIC_MATRIX_H

#include <algorithm>

#include "Types.h"
#include "Exception.h"

#ifdef __INTEL_COMPILER
#include <mkl.h>
#else

extern "C"
{
#define lapack_complex_float float _Complex // Fixes the world's most ridiculous bug
#define lapack_complex_double double _Complex
#include <lapacke.h>
#include <cblas.h>
}
#endif

template<int M, int N = 1>
class StaticMatrix
{
public:

    StaticMatrix()
    {
        std::fill(begin(), end(), 0.);
    }

    StaticMatrix(Scalar val)
    {
        std::fill(begin(), end(), val);
    }

    StaticMatrix(const std::initializer_list<Scalar> &vals)
    {
        std::copy(vals.begin(), vals.end(), vals_);
    }

    constexpr int m() const
    { return M; }

    constexpr int n() const
    { return N; }

    Scalar *data()
    { return vals_; }

    Scalar *begin()
    { return vals_; }

    Scalar *end()
    { return vals_ + M * N; }

    const Scalar *begin() const
    { return vals_; }

    const Scalar *end() const
    { return vals_ + M * N; }

    const Scalar *data() const
    { return vals_; }

    void setRow(int i, const std::initializer_list<Scalar> &vals)
    {
        std::copy(vals.begin(), vals.begin() + N, vals_ + i * N);
    }

    Scalar &operator()(int i, int j)
    { return vals_[i * N + j]; }

    Scalar operator()(int i, int j) const
    { return vals_[i * N + j]; }

    StaticMatrix<N, M> transpose() const
    {
        StaticMatrix<N, M> trans;
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                trans(j, i) = vals_[i * N + j];

        return trans;
    };

    StaticMatrix<M, 1> diag() const
    {
        StaticMatrix<M, 1> diag;
        for (int i = 0; i < M; ++i)
            diag(i, 0) = vals_[i * N + i];

        return diag;
    };

    StaticMatrix<M, N> &invert()
    {
        static_assert(M == N, "Matrix must be square.");
        lapack_int info1 = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, M, N, vals_, N, ipiv_);
        lapack_int info2 = LAPACKE_dgetri(LAPACK_ROW_MAJOR, M, vals_, N, ipiv_);

        if (info1 != 0 || info2 != 0)
            throw Exception("StaticMatrix", "invert", "matrix inversion failed.");

        return *this;
    };

    template<int K>
    void solve(StaticMatrix<M, K> &b)
    {
        static_assert(M == N, "Coefficient matrix must be square.");
        LAPACKE_dgesv(LAPACK_ROW_MAJOR, M, K, vals_, N, ipiv_, b.data(), K);
    }

    StaticMatrix<M, N> &operator*=(Scalar scalar)
    {
        for (int i = 0; i < M * N; ++i)
            vals_[i] *= scalar;

        return *this;
    }

    StaticMatrix<M, N> &operator/=(Scalar scalar)
    {
        for (int i = 0; i < M * N; ++i)
            vals_[i] /= scalar;

        return *this;
    };

private:

    Scalar vals_[M * N];
    lapack_int ipiv_[M];
};

template<int M>
StaticMatrix<M, M> eye()
{
    StaticMatrix<M, M> I;
    for (int i = 0; i < M; ++i)
        I(i, i) = 1.;
    return I;
};

template<int M>
StaticMatrix<M, M> inverse(StaticMatrix<M, M> A)
{
    A.invert();
    return A;
}

template<int M, int N>
StaticMatrix<N, M> pseudoInverse(StaticMatrix<M, N> A)
{
    StaticMatrix<M, M> I = eye<M>();
    LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', M, N, M, A.data(), N, I.data(), M);
    StaticMatrix<N, M> pInv;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            pInv(i, j) = I(i, j);

    return pInv;
}

template<int M, int K>
StaticMatrix<M, K> solve(StaticMatrix<M, M> A, StaticMatrix<M, K> b)
{
    A.solve(b);
    return b;
}

template<int M, int N>
StaticMatrix<M, N> operator-(StaticMatrix<M, N> A)
{
    std::for_each(A.begin(), A.end(), [](Scalar &a) {
        a = -a;
    });
    return A;
}

template<int M, int N>
StaticMatrix<M, N> operator*(Scalar lhs, StaticMatrix<M, N> A)
{
    A *= lhs;
    return A;
}

template<int M, int N>
StaticMatrix<M, N> operator/(StaticMatrix<M, N> A, Scalar rhs)
{
    A /= rhs;
    return A;
}

template<int M, int N, int K>
StaticMatrix<M, K> operator*(const StaticMatrix<M, N> &A, const StaticMatrix<N, K> &B)
{
    StaticMatrix<M, K> C;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                A.m(), B.n(), A.n(), 1., A.data(), A.n(),
                B.data(), B.n(), 1., C.data(), C.n());

    return C;
}

#endif
