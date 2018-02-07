#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include "TrilinosAmesosSparseMatrixSolver.h"

TrilinosAmesosSparseMatrixSolver::TrilinosAmesosSparseMatrixSolver(const Communicator &comm)
        :
        TrilinosSparseMatrixSolver(comm)
{

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

    solverName_ = parameters.get<std::string>("solver", "klu2");

    std::string filename = parameters.get<std::string>("amesosParamFile", "");

    if (!filename.empty())
        amesos2Params_ = Teuchos::getParametersFromXmlFile("case/" + filename);
    else
        amesos2Params_ = rcp(new Teuchos::ParameterList());
}