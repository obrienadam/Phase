#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include "TrilinosMueluSparseMatrixSolver.h"

TrilinosMueluSparseMatrixSolver::TrilinosMueluSparseMatrixSolver(const Communicator &comm,
                                                                 const std::string& xmlFileName)
        :
        comm_(comm)
{
    Tcomm_ = rcp(new TeuchosComm(comm.communicator()));
    mueluParameters_ = Teuchos::getParametersFromXmlFile(xmlFileName);
}

void TrilinosMueluSparseMatrixSolver::setRank(int rank)
{
    auto map = rcp(new TpetraMap(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), rank, 0, Tcomm_));
}

void TrilinosMueluSparseMatrixSolver::set(const CoefficientList& eqn)
{

}

void TrilinosMueluSparseMatrixSolver::setGuess(const Vector &x0)
{

}

void TrilinosMueluSparseMatrixSolver::setRhs(const Vector &rhs)
{

}

Scalar TrilinosMueluSparseMatrixSolver::solve()
{

}

Scalar TrilinosMueluSparseMatrixSolver::solve(const Vector &x0)
{

}

void TrilinosMueluSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{

}

void TrilinosMueluSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{

}

void TrilinosMueluSparseMatrixSolver::setMaxIters(int maxIters)
{

}

void TrilinosMueluSparseMatrixSolver::setToler(Scalar toler)
{

}

void TrilinosMueluSparseMatrixSolver::setDropToler(Scalar toler)
{

}

void TrilinosMueluSparseMatrixSolver::setFillFactor(int fill)
{

}

void TrilinosMueluSparseMatrixSolver::setMaxPreconditionerUses(int maxPreconditionerUses)
{

}

int TrilinosMueluSparseMatrixSolver::nIters() const
{

}

Scalar TrilinosMueluSparseMatrixSolver::error() const
{

}