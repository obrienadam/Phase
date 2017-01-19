#include <HYPRE_krylov.h>

#include "HypreSparseMatrixSolver.h"

HypreSparseMatrixSolver::HypreSparseMatrixSolver(const Communicator& comm)
    :
      comm_(comm)
{
    HYPRE_ParCSRBiCGSTABCreate(comm.communicator(), &solver_);
}

HypreSparseMatrixSolver::~HypreSparseMatrixSolver()
{
    HYPRE_ParCSRBiCGSTABDestroy(solver_);
}

void HypreSparseMatrixSolver::setRank(int rank)
{
    localSizes_ = comm_.allGather(rank);

    iLower_ = 0;
    for(int proc = 0; proc < comm_.rank(); ++proc)
        iLower_ += localSizes_[proc];

    iUpper_ = iLower_ + localSizes_[comm_.rank()];

    jLower_ = iLower_;
    jUpper_ = iUpper_;

    globalInds_.resize(localSizes_[comm_.rank()]);
    std::iota(globalInds_.begin(), globalInds_.end(), iLower_);

    HYPRE_IJMatrixCreate(comm_.communicator(), iLower_, iUpper_, jLower_, jUpper_, &ijMatrix_);
    HYPRE_IJMatrixSetObjectType(ijMatrix_, HYPRE_PARCSR);

    HYPRE_IJVectorCreate(comm_.communicator(), jLower_, jUpper_, &b_);
    HYPRE_IJVectorSetObjectType(b_, HYPRE_PARCSR);

    HYPRE_IJVectorCreate(comm_.communicator(), jLower_, jUpper_, &x_);
    HYPRE_IJVectorSetObjectType(x_, HYPRE_PARCSR);
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

    HYPRE_IJMatrixInitialize(ijMatrix_);
    HYPRE_IJMatrixSetValues(ijMatrix_, rows.size(), nCols.data(), rows.data(), cols.data(), vals.data());
    HYPRE_IJMatrixAssemble(ijMatrix_);
}

void HypreSparseMatrixSolver::setRhs(const Vector &rhs)
{
    std::vector<int> rows(rhs.size());

    int rowNo = jLower_;
    for(int& val: rows)
        val = rowNo++;

    HYPRE_IJVectorInitialize(b_);
    HYPRE_IJVectorSetValues(b_, rhs.size(), rows.data(), rhs.data());
    HYPRE_IJVectorAssemble(b_);

    HYPRE_IJVectorInitialize(x_);
    HYPRE_IJVectorAssemble(x_);
}

Scalar HypreSparseMatrixSolver::solve()
{
    HYPRE_IJMatrixGetObject(ijMatrix_, (void**) &csrMatrix_);
    HYPRE_IJVectorGetObject(b_, (void**) &bPar_);
    HYPRE_IJVectorGetObject(x_, (void**) &xPar_);

    HYPRE_ParCSRBiCGSTABSetup(solver_, csrMatrix_, bPar_, xPar_);
    HYPRE_ParCSRBiCGSTABSolve(solver_, csrMatrix_, bPar_, xPar_);

    HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver_, &toler_);
    HYPRE_BiCGSTABGetNumIterations(solver_, &nIters_);

    return toler_;
}

void HypreSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{
    std::vector<int> inds(field.grid.activeCells().size());
    std::vector<Scalar> vals(field.grid.activeCells().size());

//    std::transform(field.grid.activeCells().begin(), field.grid.activeCells().end(),
//                   inds.begin(), [](const Cell& cell)->int{ return cell.globalIndex(); });

    HYPRE_IJVectorGetValues(x_, inds.size(), inds.data(), vals.data());

    auto begin = vals.begin();
    for(const Cell& cell: field.grid.activeCells())
        field(cell) = *(begin++);
}

void HypreSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{
    std::vector<int> inds(2*field.grid.activeCells().size());
    std::vector<Scalar> vals(2*field.grid.activeCells().size());

    Size nActiveCells = field.grid.activeCells().size();

//    std::transform(field.grid.activeCells().begin(), field.grid.activeCells().end(),
//                   inds.begin(), [](const Cell& cell)->int{ return cell.globalIndex(); });

//    std::transform(field.grid.activeCells().begin(), field.grid.activeCells().end(),
//                   inds.begin() + nActiveCells,
//                   [nActiveCells](const Cell& cell)->int{ return cell.globalIndex() + nActiveCells; });

    HYPRE_IJVectorGetValues(x_, inds.size(), inds.data(), vals.data());

    auto begin = vals.begin();
    for(const Cell& cell: field.grid.activeCells())
    {
        field(cell).x = *begin;
        field(cell).y = *(begin + nActiveCells);
        ++begin;
    }
}
