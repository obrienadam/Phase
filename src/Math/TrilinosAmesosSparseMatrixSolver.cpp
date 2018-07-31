#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <TpetraExt_MatrixMatrix.hpp>

#include "TrilinosAmesosSparseMatrixSolver.h"

TrilinosAmesosSparseMatrixSolver::TrilinosAmesosSparseMatrixSolver(const Communicator &comm, const std::string &solverName)
    :
      TrilinosSparseMatrixSolver(comm)
{
    solverName_ = solverName;
    amesos2Params_ = rcp(new Teuchos::ParameterList());
}

TrilinosAmesosSparseMatrixSolver::TrilinosAmesosSparseMatrixSolver(const Communicator &comm, Tpetra::ProfileType pftype, const std::string &solverName)
    :
      TrilinosSparseMatrixSolver(comm, pftype)
{
    solverName_ = solverName;
    amesos2Params_ = rcp(new Teuchos::ParameterList());
}

Scalar TrilinosAmesosSparseMatrixSolver::solve()
{
    if (solver_.is_null())
    {
        solver_ = Amesos2::create<TpetraCrsMatrix, TpetraMultiVector>(solverName_, mat_, x_, b_);
        solver_->setParameters(amesos2Params_);
    }
    else
    {
        solver_->setA(mat_);
        solver_->setX(x_);
        solver_->setB(b_);
    }

    solver_->symbolicFactorization().numericFactorization().solve();

    return error();
}

void TrilinosAmesosSparseMatrixSolver::setup(const boost::property_tree::ptree &parameters)
{
    using namespace Teuchos;

    solverName_ = parameters.get<std::string>("solver", solverName_);

    std::string filename = parameters.get<std::string>("amesosParamFile", "");

    if (!filename.empty())
        amesos2Params_ = Teuchos::getParametersFromXmlFile("case/" + filename);
}

Matrix TrilinosAmesosSparseMatrixSolver::pseudoInverse() const
{
    TpetraMultiVector I(mat_->getRowMap(), mat_->getGlobalNumRows(), true); //- m x m identity matrix

    Size nRowsLocal = I.getLocalLength();

    auto view = I.get2dViewNonConst();

    for(int i = 0; i < nRowsLocal; ++i)
        view[i][i] = 1.;

    auto A = Teuchos::rcp(new TpetraCrsMatrix(mat_->getDomainMap(), 12)); //- matrix will be n x n
    auto x = Teuchos::rcp(new TpetraMultiVector(mat_->getDomainMap(), mat_->getGlobalNumRows(), false)); //- soln will be n x m
    auto b = Teuchos::rcp(new TpetraMultiVector(mat_->getDomainMap(), mat_->getGlobalNumRows(), false)); //- rhs will be n x m

    Tpetra::MatrixMatrix::Multiply(*mat_, true, *mat_, false, *A, true);
    mat_->apply(I, *b, Teuchos::TRANS);

    auto solver = Amesos2::create<TpetraCrsMatrix, TpetraMultiVector>("klu2", A, x, b);

    solver->symbolicFactorization().numericFactorization().solve();

    nRowsLocal = x->getLocalLength();
    Size nCols = x->getNumVectors();
    view = x->get2dViewNonConst();

    Matrix soln(nRowsLocal, nCols);

    for(int i = 0; i < nRowsLocal; ++i)
        for(int j = 0; j < nCols; ++j)
            soln(i, j) = view[i][j];

    return soln;
}
