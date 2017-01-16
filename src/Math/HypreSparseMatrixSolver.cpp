#include <HYPRE_krylov.h>

#include "HypreSparseMatrixSolver.h"

HypreSparseMatrixSolver::HypreSparseMatrixSolver(const Communicator& comm)
    :
      comm_(comm)
{

}

void HypreSparseMatrixSolver::setRank(int rank)
{
    localSizes_ = comm_.allGather(rank);

    iLower_ = 0;
    for(int proc = 0; proc < comm_.rank(); ++proc)
        iLower_ += localSizes_[proc];

    iUpper_ = iLower_ + localSizes_[comm_.rank()];

    jLower_ = iLower_;
    jUpper_ = jUpper_;

    HYPRE_IJMatrixDestroy(ijMatrix_);
    HYPRE_IJMatrixCreate(comm_.communicator(), iLower_, iUpper_, jLower_, jUpper_, &ijMatrix_);
    HYPRE_IJMatrixSetObjectType(ijMatrix_, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(ijMatrix_);

    HYPRE_IJVectorDestroy(b_);
    HYPRE_IJVectorCreate(comm_.communicator(), jLower_, jUpper_, &b_);
    HYPRE_IJVectorInitialize(b_);

    HYPRE_IJVectorDestroy(x_);
    HYPRE_IJVectorCreate(comm_.communicator(), jLower_, jUpper_, &x_);
    HYPRE_IJVectorInitialize(x_);
}

void HypreSparseMatrixSolver::set(const SparseMatrixSolver::CoefficientList &eqn)
{
    std::vector<int> nCols;
    std::vector<int> rows;
    nCols.reserve(eqn.size());
    rows.reserve(eqn.size());

    Size nEntries = 0;
    int rowNo = iLower_;
    for(const auto& row: eqn)
    {
        nCols.push_back(row.size());
        nEntries += row.size();
        rows.push_back(rowNo++);
    }

    std::vector<int> cols;
    std::vector<Scalar> vals;
    cols.reserve(nEntries);
    vals.reserve(nEntries);

    for(const auto& row: eqn)
        for(const auto& entry: row)
        {
            cols.push_back(entry.first);
            vals.push_back(entry.second);
        }

    HYPRE_IJMatrixSetValues(ijMatrix_, rows.size(), nCols.data(), rows.data(), cols.data(), vals.data());
    HYPRE_IJMatrixAssemble(ijMatrix_);
}

void HypreSparseMatrixSolver::setRhs(const Vector &rhs)
{
    std::vector<int> rows(rhs.size());

    int rowNo = jLower_;
    for(int& val: rows)
        val = rowNo++;

    HYPRE_IJVectorSetValues(b_, rhs.size(), rows.data(), rhs.data());
    HYPRE_IJVectorAssemble(b_);
    HYPRE_IJVectorAssemble(x_);
}

Scalar HypreSparseMatrixSolver::solve()
{
    HYPRE_IJMatrixGetObject(ijMatrix_, (void**) &csrMatrix_);
    HYPRE_IJVectorGetObject(b_, (void**) &bPar_);
    HYPRE_IJVectorGetObject(x_, (void**) &xPar_);

    HYPRE_BiCGSTABSetup(solver_, (HYPRE_Matrix)csrMatrix_, (HYPRE_Vector)bPar_, (HYPRE_Vector)xPar_);
    HYPRE_BiCGSTABSolve(solver_, (HYPRE_Matrix)csrMatrix_, (HYPRE_Vector)bPar_, (HYPRE_Vector)xPar_);

    HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver_, &toler_);
    HYPRE_BiCGSTABGetNumIterations(solver_, &nIters_);

    HYPRE_BiCGSTABDestroy(solver_);

    return toler_;
}

void HypreSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{

}

void HypreSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{

}
