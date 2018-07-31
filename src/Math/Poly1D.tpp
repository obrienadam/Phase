#include "Poly1D.h"
#include "StaticMatrix.h"

template<Size N>
Poly1D<N>::Poly1D(const std::initializer_list<Scalar> &x,
                  const std::initializer_list<Scalar> &rhs)
{
    StaticMatrix<N, N> A;

    auto itr = x.begin();
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            A(i, j) = std::pow(*(itr++), i);

    _coeffs = rhs;
    A.solve(_coeffs);
}

template<Size N>
Scalar Poly1D<N>::operator()(Scalar x) const
{
    StaticMatrix<1, N> v;
    for(int i = 0; i < N; ++i)
        v(0, i) = std::pow(x, i);

    return v * _coeffs;
}
