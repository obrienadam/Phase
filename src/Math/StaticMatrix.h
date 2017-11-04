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
        std::fill(vals_, vals_ + M * N, 0.);
    }

    StaticMatrix(const std::initializer_list<Scalar> &vals)
    {
        std::transform(vals.begin(), vals.end(), vals_, [](Scalar val) { return val; });
    }

    constexpr int m() const
    { return M; }

    constexpr int n() const
    { return N; }

    Scalar *data()
    { return vals_; }

    const Scalar *begin() const
    { return vals_; }

    const Scalar *end() const
    { return vals_ + M * N; }

    const Scalar *data() const
    { return vals_; }

    Scalar &operator()(int i, int j)
    { return vals_[i * N + j]; }

    Scalar operator()(int i, int j) const
    { return vals_[i * N + j]; }

    StaticMatrix<N, M> transpose() const
    {
        StaticMatrix<N, M> matT;
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                matT(j, i) = vals_[i * N + j];

        return matT;
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

        for (int i = 0; i < M; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                std::cout << vals_[i * N + j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";

        lapack_int info1 = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, M, N, vals_, N, ipiv_);
        lapack_int info2 = LAPACKE_dgetri(LAPACK_ROW_MAJOR, M, vals_, N, ipiv_);

        if (info1 != 0 || info2 != 0)
            throw Exception("StaticMatrix", "invert", "inversion failed, matrix is singular to working precision.");

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
    int ipiv_[M];
};

template<int M, int N>
StaticMatrix<M, N> inverse(StaticMatrix<M, N> A)
{
    A.invert();
    return A;
}

template<int M, int N>
StaticMatrix<M, M> pseudoInverse(StaticMatrix<M, N> A)
{
    return inverse(A * A.transpose()) * A;
}

template<int M, int N, int K>
StaticMatrix<M, K> solve(StaticMatrix<M, N> A, StaticMatrix<M, K> b)
{
    A.solve(b);
    return b;
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
