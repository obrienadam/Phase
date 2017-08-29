#include <stdio.h>

#include "Equation.h"
#include "Exception.h"
#include "SparseMatrixSolver.h"

template<class T>
Equation<T>::Equation(const Input &input,
                      FiniteVolumeField<T> &field,
                      const std::string &name)
        :
        Equation<T>::Equation(field, name)
{
    configureSparseSolver(input, field.grid().comm());
}

template<class T>
void Equation<T>::clear()
{
    for (auto &row: coeffs_)
        row.clear();

    sources_.zero();
}

template<class T>
Equation<T> &Equation<T>::operator=(const Equation<T> &rhs)
{
    if (this == &rhs)
        return *this;
    else if (&field_.grid() != &rhs.field_.grid())
        throw Exception("Equation<T>", "operator=", "cannot copy equations defined for different fields.");

    nLocalActiveCells_ = rhs.nLocalActiveCells_;
    nGlobalActiveCells_ = rhs.nGlobalActiveCells_;
    coeffs_ = rhs.coeffs_;
    sources_ = rhs.sources_;

    if(rhs.indexMap_)
        indexMap_ = rhs.indexMap_;

    if (rhs.spSolver_) // Prevent a sparse solver from being accidently destroyed if the rhs solver doesn't exist
        spSolver_ = rhs.spSolver_;

    return *this;
}

template<class T>
Equation<T> &Equation<T>::operator=(Equation<T> &&rhs)
{
    if (&field_ != &rhs.field_)
        throw Exception("Equation<T>", "operator=", "cannot copy equations defined for different fields.");

    nLocalActiveCells_ = rhs.nLocalActiveCells_;
    nGlobalActiveCells_ = rhs.nGlobalActiveCells_;
    coeffs_ = std::move(rhs.coeffs_);
    sources_ = std::move(rhs.sources_);

    if (rhs.spSolver_)
        spSolver_ = rhs.spSolver_;

    return *this;
}

template<class T>
Equation<T> &Equation<T>::operator+=(const Equation<T> &rhs)
{
    for (int i = 0; i < rhs.coeffs_.size(); ++i)
        for (const auto &entry: rhs.coeffs_[i])
            addValue(i, entry.first, entry.second);

    sources_ += rhs.sources_;

    return *this;
}

template<class T>
Equation<T> &Equation<T>::operator-=(const Equation<T> &rhs)
{
    for (int i = 0; i < rhs.coeffs_.size(); ++i)
        for (const auto &entry: rhs.coeffs_[i])
            add(i, entry.first, -entry.second);

    sources_ -= rhs.sources_;

    return *this;
}

template<class T>
Equation<T> &Equation<T>::operator*=(Scalar rhs)
{
    for (int i = 0; i < coeffs_.size(); ++i)
        for (const auto &entry: coeffs_[i])
            set(i, entry.first, rhs * entry.second);

    sources_ *= rhs;

    return *this;
}

template<class T>
Equation<T> &Equation<T>::operator/=(const ScalarFiniteVolumeField &rhs)
{
    for(const Cell& cell: rhs.grid().localActiveCells())
    {
        Index i = cell.index(0);
        Scalar val = rhs(cell);

        for(auto& coeff: coeffs_[i])
            coeff.second /= val;

        sources_[i] /= val;
    }

    return *this;
}

template<class T>
Equation<T> &Equation<T>::operator==(Scalar rhs)
{
    if(rhs == 0.) // for efficiency
        return *this;

    for (Scalar &src: sources_)
        src -= rhs;

    return *this;
}

template<class T>
Equation<T> &Equation<T>::operator==(const Equation<T> &rhs)
{
    for (int i = 0; i < coeffs_.size(); ++i)
        for (const auto &entry: rhs.coeffs_[i])
            addValue(i, entry.first, -entry.second);

    sources_ -= rhs.sources_;

    return *this;
}

template<class T>
Equation<T> &Equation<T>::operator==(const FiniteVolumeField<T> &rhs)
{
    for (const Cell &cell: rhs.grid().localActiveCells())
        sources_(cell.index(0)) -= rhs(cell);

    return *this;
}

template<class T>
void Equation<T>::setSparseSolver(std::shared_ptr<SparseMatrixSolver> &spSolver)
{
    spSolver_ = spSolver;
}

template<class T>
void Equation<T>::configureSparseSolver(const Input &input, const Communicator &comm)
{
    std::string lib = input.caseInput().get<std::string>("LinearAlgebra." + name + ".lib", "Eigen3");
    boost::algorithm::to_lower(lib);

    if (lib == "eigen" || lib == "eigen3")
        spSolver_ = std::make_shared<EigenSparseMatrixSolver>();
    else if(lib == "trilinos" || lib == "belos")
        spSolver_ = std::make_shared<TrilinosBelosSparseMatrixSolver>(comm);
    //else if(lib == "muelu")
    //    spSolver_ = std::make_shared<TrilinosMueluSparseMatrixSolver>(comm);
    else
        throw Exception("Equation<T>", "configureSparseSolver", "unrecognized sparse solver lib \"" + lib + "\".");

    if (comm.nProcs() > 1 && !spSolver_->supportsMPI())
        throw Exception("Equation<T>", "configureSparseSolver", "equation \"" + name + "\", lib \"" + lib +
                                                                "\" does not support multiple processes in its current configuration.");

    spSolver_->setup(input.caseInput().get_child("LinearAlgebra." + name));

    comm.printf("Initialized sparse matrix solver for equation \"%s\" using lib%s.\n", name.c_str(), lib.c_str());
}

template<class T>
Scalar Equation<T>::solve()
{
    if (!spSolver_)
        throw Exception("Equation<T>", "solve",
                        "must allocate a SparseMatrixSolver object before attempting to solve.");

    nLocalActiveCells_ = field_.grid().nLocalActiveCells();
    nGlobalActiveCells_ = field_.grid().nActiveCellsGlobal();

    spSolver_->setRank(getRank());
    spSolver_->set(coeffs_);
    spSolver_->setRhs(-sources_);
    spSolver_->solve();
    spSolver_->mapSolution(field_);

    spSolver_->printStatus("Equation " + name + ":");

    return spSolver_->error();
}

//- Private methods

template<class T>
void Equation<T>::setValue(Index i, Index j, Scalar val)
{
    for (auto &entry: coeffs_[i])
    {
        if (entry.first == j)
        {
            entry.second = val;
            return;
        }
    }

    coeffs_[i].push_back(std::make_pair(j, val));
}

template<class T>
void Equation<T>::addValue(Index i, Index j, Scalar val)
{
    for (auto &entry: coeffs_[i])
    {
        if (entry.first == j)
        {
            entry.second += val;
            return;
        }
    }

    coeffs_[i].push_back(std::make_pair(j, val));
}

template<class T>
Scalar &Equation<T>::coeffRef(Index i, Index j)
{
    for (auto &entry: coeffs_[i])
    {
        if (entry.first == j)
            return entry.second;
    }

    throw Exception("Equation<T>", "coeffRef", "requested coefficient does not exist.");
}

//- External functions

template<class T>
Equation<T> operator+(Equation<T> lhs, const Equation<T> &rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
Equation<T> operator-(Equation<T> lhs, const Equation<T> &rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class T>
Equation<T> operator+(Equation<T> lhs, const FiniteVolumeField<T> &rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
Equation<T> operator+(const FiniteVolumeField<T> &lhs, Equation<T> rhs)
{
    rhs += lhs;
    return rhs;
}

template<class T>
Equation<T> operator-(Equation<T> lhs, const FiniteVolumeField<T> &rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class T>
Equation<T> operator-(const FiniteVolumeField<T> &lhs, Equation<T> rhs)
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
