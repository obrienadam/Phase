#include <stdio.h>

#include "Equation.h"
#include "Exception.h"

template<class T>
Equation<T>::Equation(const Input& input, T& field, const std::string& name)
    :
      Equation<T>::Equation(field, name)
{
    solver_.setMaxIterations(input.caseInput().get<int>("LinearAlgebra." + name + ".maxIterations", 500));
    solver_.setTolerance(input.caseInput().get<Scalar>("LinearAlgebra." + name + ".tolerance", 1e-10));
    solver_.preconditioner().setFillfactor(input.caseInput().get<int>("LinearAlgebra." + name + ".iluFill", 3));
}

template<class T>
Equation<T>::Equation(const Equation<T>& other)
    :
      coeffs_(other.coeffs_),
      boundaries_(other.boundaries_),
      sources_(other.sources_),
      field_(other.field_)
{

}

template<class T>
void Equation<T>::set(int i, int j, Scalar val)
{
    for(auto& coeff: coeffs_[i])
    {
        if(coeff.first == j)
        {
            coeff.second = val;
            return;
        }
    }

    coeffs_[i].push_back(std::make_pair(j, val));
}

template <class T>
void Equation<T>::add(int i, int j, Scalar val)
{
    for(auto& coeff: coeffs_[i])
    {
        if(coeff.first == j)
        {
            coeff.second += val;
            return;
        }
    }

    coeffs_[i].push_back(std::make_pair(j, val));
}

template<class T>
Scalar& Equation<T>::getRef(int i, int j)
{
    for(auto& coeff: coeffs_[i])
    {
        if(coeff.first == j)
            return coeff.second;
    }

    throw Exception("Equation<T>", "getRef", "no such entry.");
}

template<class T>
Scalar Equation<T>::get(int i, int j) const
{
    for(const auto& coeff: coeffs_[i])
    {
        if(coeff.first == j)
            return coeff.second;
    }

    return 0.;
}

template<class T>
void Equation<T>::clear()
{
    for(auto& coeff: coeffs_)
        coeff.clear();
}

template<class T>
Scalar Equation<T>::solve()
{
    std::vector<Triplet> triplets;
    for(int i = 0; i < coeffs_.size(); ++i)
        for(const auto& entry: coeffs_[i])
            triplets.push_back(Triplet(i, entry.first, entry.second));

    spMat_.setFromTriplets(triplets.begin(), triplets.end());

    solver_.compute(spMat_);
    field_ = solver_.solve(boundaries_ + sources_);

    printf("Solved %s equation, error = %lf, number of iterations = %d\n", name.c_str(), error(), iterations());
    return error();
}

template<class T>
Scalar Equation<T>::solve(const SparseVector &x0)
{
    std::vector<Triplet> triplets;
    for(int i = 0; i < coeffs_.size(); ++i)
        for(const auto& entry: coeffs_[i])
            triplets.push_back(Triplet(i, entry.first, entry.second));

    spMat_.setFromTriplets(triplets.begin(), triplets.end());

    solver_.compute(spMat_);
    field_ = solver_.solveWithGuess(boundaries_ + sources_, x0);

    printf("Solved %s equation, error = %lf, number of iterations = %d\n", name.c_str(), error(), iterations());
    return error();
}

template<class T>
Equation<T>& Equation<T>::operator +=(const Equation<T>& rhs)
{
    for(int i = 0; i < rhs.coeffs_.size(); ++i)
        for(const auto& entry: rhs.coeffs_[i])
            add(i, entry.first, entry.second);

    boundaries_ += rhs.boundaries_;
    sources_ += rhs.sources_;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator -=(const Equation<T>& rhs)
{
    for(int i = 0; i < rhs.coeffs_.size(); ++i)
        for(const auto& entry: rhs.coeffs_[i])
            add(i, entry.first, -entry.second);

    boundaries_ -= rhs.boundaries_;
    sources_ -= rhs.sources_;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator =(const Equation<T>& rhs)
{
    if(this != &rhs)
    {
        coeffs_ = rhs.coeffs_;
        boundaries_ = rhs.boundaries_;
        sources_ = rhs.sources_;
    }

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator *=(Scalar rhs)
{
    for(int i = 0; i < coeffs_.size(); ++i)
        for(const auto& entry: coeffs_[i])
            set(i, entry.first, rhs*entry.second);

    boundaries_ *= rhs;
    sources_ *= rhs;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator==(Scalar rhs)
{
    for(int i = 0, end = boundaries_.rows(); i < end; ++i)
        sources_(i) += rhs;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator==(const Equation<T>& rhs)
{
    for(int i = 0; i < rhs.coeffs_.size(); ++i)
        for(const auto& entry: rhs.coeffs_[i])
            add(i, entry.first, -entry.second);

    boundaries_ -= rhs.boundaries_;
    sources_ = rhs.sources_ - sources_;
    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator ==(const T& rhs)
{
    for(const Cell& cell: rhs.grid.fluidCells())
        sources_(cell.globalIndex()) += rhs[cell.id()];

    return *this;
}

//- External functions

template<class T>
Equation<T> operator+(Equation<T> lhs, const Equation<T>& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
Equation<T> operator-(Equation<T> lhs, const Equation<T>& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class T>
Equation<T> operator+(Equation<T> lhs, const T& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
Equation<T> operator+(const T& lhs, Equation<T> rhs)
{
    rhs += lhs;
    return rhs;
}

template<class T>
Equation<T> operator-(Equation<T> lhs, const T& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class T>
Equation<T> operator-(const T& lhs, Equation<T> rhs)
{
    rhs -= lhs;
    return rhs;
}

template<class T>
Equation<T> operator*(Equation<T> lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

template<class T>
Equation<T> operator*(Scalar lhs, Equation<T> rhs)
{
    rhs *= lhs;
    return rhs;
}
