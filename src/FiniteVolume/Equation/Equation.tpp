#include <stdio.h>

#include "Equation.h"
#include "Exception.h"
#include "SparseMatrixSolver.h"

template<class T>
Equation<T>::Equation(const Input& input,
                      const Communicator& comm,
                      FiniteVolumeField<T> &field,
                      const std::string& name)
    :
      Equation<T>::Equation(field, name)
{
    configureSparseSolver(input, comm);
}

template<class T>
void Equation<T>::clear()
{
    for(auto& coeff: coeffs_)
        coeff.clear();
}

template<class T>
Equation<T> &Equation<T>::operator=(const Equation<T> &rhs)
{
    if(this == &rhs)
        return *this;
    else if(&field_ != &rhs.field_)
        throw Exception("Equation<T>", "operator=", "cannot copy equations defined for different fields.");

    coeffs_ = rhs.coeffs_;
    boundaries_ = rhs.boundaries_;
    sources_ = rhs.sources_;

    if(rhs.spSolver_) // Prevent a sparse solver from being accidently destroyed if the rhs solver doesn't exist
        spSolver_ = rhs.spSolver_;

    return *this;
}

template<class T>
Equation<T> &Equation<T>::operator=(Equation<T> &&rhs)
{
    if(&field_ != &rhs.field_)
        throw Exception("Equation<T>", "operator=", "cannot copy equations defined for different fields.");

    coeffs_ = std::move(rhs.coeffs_);
    boundaries_ = std::move(rhs.boundaries_);
    sources_ = std::move(sources_);

    if(rhs.spSolver_)
        spSolver_ = rhs.spSolver_;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator +=(const Equation<T>& rhs)
{
    for(int i = 0; i < rhs.coeffs_.size(); ++i)
        for(const auto& entry: rhs.coeffs_[i])
            addValue(i, entry.first, entry.second);

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
            addValue(i, entry.first, -entry.second);

    boundaries_ -= rhs.boundaries_;
    sources_ = rhs.sources_ - sources_;
    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator ==(const FiniteVolumeField<T>& rhs)
{
    for(const Cell& cell: rhs.grid.localActiveCells())
        sources_(cell.localIndex()) += rhs(cell);

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

    if(lib == "eigen" || lib == "eigen3")
#ifdef PHASE_USE_EIGEN
        spSolver_ = std::shared_ptr<EigenSparseMatrixSolver>(new EigenSparseMatrixSolver());
#else
        throw Exception("Equation<T>", "configureSparseMatrixSolver", "Phase was built without Eigen support. Recompile with \"-DPHASE_USE_EIGEN=ON\".");
#endif
    else if(lib == "hypre")
#ifdef PHASE_USE_HYPRE
        spSolver_ = std::shared_ptr<HypreSparseMatrixSolver>(new HypreSparseMatrixSolver(comm));
#else
        throw Exception("Equation<T>", "configureSparseMatrixSolver", "Phase was built without HYPRE support. Recompile with \"-DPHASE_USE_HYPRE=ON\".");
#endif
    else if(lib == "petsc")
    {
#ifdef PHASE_USE_PETSC
        std::string precon = input.caseInput().get<std::string>("LinearAlgebra." + name + ".preconditioner", "sor");
        spSolver_ = std::shared_ptr<PetscSparseMatrixSolver>(new PetscSparseMatrixSolver(comm, precon));
#else
        throw Exception("Equation<T>", "configureSparseMatrixSolver", "Phase was built without Petsc support. Recompile with \"-DPHASE_USE_PETSC=ON\".");
#endif
    }
    else
        throw Exception("Equation<T>", "configureSparseSolver", "unrecognized sparse solver lib \"" + lib + "\".");

    if(comm.nProcs() > 1 && !spSolver_->supportsMPI())
        throw Exception("Equation<T>", "configureSparseSolver", "equation \"" + name + "\", lib \"" + lib + "\" does not support multiple processes in its current configuration.");

    spSolver_->setMaxIters(input.caseInput().get<int>("LinearAlgebra." + name + ".maxIterations", 500));
    spSolver_->setToler(input.caseInput().get<Scalar>("LinearAlgebra." + name + ".tolerance", 1e-6));
    spSolver_->setFillFactor(input.caseInput().get<int>("LinearAlgebra." + name + ".iluFill", 2));
    spSolver_->setDropToler(input.caseInput().get<Scalar>("LinearAlgebra." + name + ".dropTolerance", 0));
    spSolver_->setMaxPreconditionerUses(input.caseInput().get<int>("LinearAlgebra." + name + ".maxPreconditionerUses", 1));

    comm.printf("Initialized sparse matrix solver for equation \"%s\" using lib%s.\n", name.c_str(), lib.c_str());
}

template<class T>
Scalar Equation<T>::solve()
{
    if(!spSolver_)
        throw Exception("Equation<T>", "solve", "must allocate a SparseMatrixSolver object before attempting to solve.");

    nActiveCells_ = field_.grid.nLocalActiveCells();

    spSolver_->setRank(getRank());
    spSolver_->set(coeffs_);
    spSolver_->setRhs(boundaries_ + sources_);
    spSolver_->solve();
    spSolver_->mapSolution(field_);

    spSolver_->printStatus("Equation " + name + ":");

    return spSolver_->error();
}

template<class T>
Scalar Equation<T>::solveWithGuess()
{
    if(!spSolver_)
        throw Exception("Equation<T>", "solve", "must allocate a SparseMatrixSolver object before attempting to solve.");

    nActiveCells_ = field_.grid.nLocalActiveCells();

    spSolver_->setRank(getRank());
    spSolver_->set(coeffs_);
    spSolver_->setRhs(boundaries_ + sources_);
    spSolver_->solve(field_.vectorize());
    spSolver_->mapSolution(field_);

    spSolver_->printStatus("Equation " + name + ":");

    return spSolver_->error();
}

//- Private methods

template<class T>
void Equation<T>::setValue(Index i, Index j, Scalar val)
{
    for(auto& entry: coeffs_[i])
    {
        if(entry.first == j)
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
    for(auto& entry: coeffs_[i])
    {
        if(entry.first == j)
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
    for(auto& entry: coeffs_[i])
    {
        if(entry.first == j)
            return entry.second;
    }

    throw Exception("Equation<T>", "coeff", "requested coefficient does not exist.");
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
Equation<T> operator+(Equation<T> lhs, const FiniteVolumeField<T>& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
Equation<T> operator+(const FiniteVolumeField<T>& lhs, Equation<T> rhs)
{
    rhs += lhs;
    return rhs;
}

template<class T>
Equation<T> operator-(Equation<T> lhs, const FiniteVolumeField<T>& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class T>
Equation<T> operator-(const FiniteVolumeField<T>& lhs, Equation<T> rhs)
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
