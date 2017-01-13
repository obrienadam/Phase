#include <stdio.h>

#include "Equation.h"
#include "Exception.h"

template<class T>
Equation<T>::Equation(const Input& input, T& field, const std::string& name, std::shared_ptr<SparseMatrixSolver> &&spSolver)
    :
      Equation<T>::Equation(field, name)
{
    spSolver_ = spSolver;
    configureSparseSolver(input);
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
    for(int i = 0, end = boundaries_.size(); i < end; ++i)
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

template<class T>
void Equation<T>::setSparseSolver(std::shared_ptr<SparseMatrixSolver> &spSolver)
{
    spSolver_ = spSolver;
}

template<class T>
void Equation<T>::configureSparseSolver(const Input &input)
{
    if(!spSolver_)
        throw Exception("Equation<T>", "configureSparseSolver", "must allocate a SparseMatrixSolver object before attempting to configure.");

    spSolver_->setMaxIters(input.caseInput().get<int>("LinearAlgebra." + name + ".maxIterations", 500));
    spSolver_->setToler(input.caseInput().get<Scalar>("LinearAlgebra." + name + ".tolerance", 1e-10));
    spSolver_->setFillFactor(input.caseInput().get<int>("LinearAlgebra." + name + ".iluFill", 3));
}

template<class T>
Scalar Equation<T>::solve()
{
    if(!spSolver_)
        throw Exception("Equation<T>", "solve", "must allocate a SparseMatrixSolver object before attempting to solve.");

    spSolver_->setRank(field_.dimension()*field_.grid.nActiveCells());
    spSolver_->set(coeffs_);
    spSolver_->setRhs(boundaries_ + sources_);
    spSolver_->solve();
    spSolver_->mapSolution(field_);

    printf("Solved equation \"%s\" in %d iterations, error = %lf\n", name.c_str(), spSolver_->nIters(), spSolver_->error());
    return spSolver_->error();
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
