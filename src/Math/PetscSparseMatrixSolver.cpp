#include "PetscSparseMatrixSolver.h"

int PetscSparseMatrixSolver::solversActive_ = 0;

PetscSparseMatrixSolver::PetscSparseMatrixSolver(const Communicator &comm, const std::string &preconditioner)
    :
      comm_(comm)
{
    if(!PetscInitializeCalled)
        PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

    error_ = MatCreate(comm_.communicator(), &A_);
    error_ = MatSetType(A_, MATAIJ);
    error_ = VecCreate(comm_.communicator(), &b_);
    error_ = VecSetType(b_, VECMPI);
    error_ = VecCreate(comm_.communicator(), &x_);
    error_ = VecSetType(x_, VECMPI);

    KSPCreate(comm_.communicator(), &solver_);
    KSPSetOperators(solver_, A_, A_);
    KSPSetType(solver_, KSPBCGS);

    setPreconditioner(preconditioner);
    PetscViewerCreate(comm.communicator(), &viewer_);

    ++solversActive_;
}

PetscSparseMatrixSolver::~PetscSparseMatrixSolver()
{
    MatDestroy(&A_);
    VecDestroy(&x_);
    VecDestroy(&b_);
    KSPDestroy(&solver_);
    PetscViewerDestroy(&viewer_);
}

void PetscSparseMatrixSolver::setRank(int rank)
{
    if(rank != iUpper_ - iLower_)
    {
        error_ = MatSetSizes(A_, rank, rank, PETSC_DETERMINE, PETSC_DETERMINE);
        error_ = MatSeqAIJSetPreallocation(A_, 5, PETSC_NULL);
        error_ = MatMPIAIJSetPreallocation(A_, 5, PETSC_NULL, 2, PETSC_NULL);
        error_ = MatGetOwnershipRange(A_, &iLower_, &iUpper_);
        error_ = VecSetSizes(x_, rank, PETSC_DECIDE);
        error_ = VecSetSizes(b_, rank, PETSC_DECIDE);
    }
}

void PetscSparseMatrixSolver::set(const SparseMatrixSolver::CoefficientList &eqn)
{
    PetscInt rowNo = iLower_;
    for(const Row& row: eqn)
    {
        for(const Entry& entry: row)
            MatSetValue(A_, rowNo, entry.first, entry.second, INSERT_VALUES);

        rowNo++;
    }

    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
}

void PetscSparseMatrixSolver::setRhs(const Vector &rhs)
{
    PetscInt rowNo = iLower_;
    for(Scalar val: rhs)
        VecSetValue(b_, rowNo++, val, INSERT_VALUES);

    VecAssemblyBegin(b_);
    VecAssemblyEnd(b_);
}

Scalar PetscSparseMatrixSolver::solve()
{
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);

    KSPDestroy(&solver_);
    KSPCreate(comm_.communicator(), &solver_); // Unfortunately necessary to fix memory leak!
    KSPSetOperators(solver_, A_, A_);
    KSPSetType(solver_, KSPBCGS);

    setPreconditioner("pilut");

    KSPSolve(solver_, b_, x_);
    return error();
}

void PetscSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{
    PetscInt row = iLower_;
    for(const Cell& cell: field.grid.activeCells())
    {
        VecGetValues(x_, 1, &row, &(field(cell)));
        ++row;
    }
}

void PetscSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{
    PetscInt rowX = iLower_;
    PetscInt rowY = rowX + field.grid.nActiveCells();
    for(const Cell& cell: field.grid.activeCells())
    {
        VecGetValues(x_, 1, &rowX, &(field(cell).x));
        VecGetValues(x_, 1, &rowY, &(field(cell).y));

        ++rowX;
        ++rowY;
    }
}

void PetscSparseMatrixSolver::setPreconditioner(const std::string &preconditioner)
{
    KSPGetPC(solver_, &precon_);

    if(preconditioner == "pilut")
    {
        PCSetType(precon_, PCHYPRE);
        PCHYPRESetType(precon_, "pilut");
        PetscOptionsSetValue(NULL, "-pc_hypre_pilut_factorrowsize","100");
        PetscOptionsSetValue(NULL, "-pc_hypre_pilut_tol", "0.001");
        PCSetFromOptions(precon_);
    }
    else if(preconditioner == "boomeramg")
    {
        PCSetType(precon_, PCHYPRE);
        PCHYPRESetType(precon_, "boomeramg");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_down", "sequential-Gauss-Seidel");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_up", "sequential-Gauss-Seidel");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_coarse", "Gaussian-elimination");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_max_iter", "5");
        PCSetFromOptions(precon_);
    }
    else if(preconditioner == "sor")
    {
        PCSetType(precon_, PCSOR);
    }
    else if(preconditioner == "lu")
    {
        PCSetType(precon_, PCLU);
        KSPSetType(solver_, KSPPREONLY);
    }
    else
        throw Exception("PetscSparseMatrixSolver", "setPreconditioner", "unrecognized preconditioner type \"" + preconditioner + "\".");
}

void PetscSparseMatrixSolver::setMaxIters(int maxIters)
{
    maxIters_ = maxIters;
    KSPSetTolerances(solver_, rTol_, absTol_, PETSC_DEFAULT, maxIters_);
}

void PetscSparseMatrixSolver::setToler(Scalar toler)
{
    rTol_ = absTol_ = toler;
    KSPSetTolerances(solver_, rTol_, absTol_, PETSC_DEFAULT, maxIters_);
}

void PetscSparseMatrixSolver::setDropToler(Scalar toler)
{
    //PetscOptionsSetValue(NULL, "-pc_hypre_pilut_tol", std::to_string(toler).c_str());
    //PCSetFromOptions(precon_);
}

void PetscSparseMatrixSolver::setFillFactor(int fill)
{
    //PetscOptionsSetValue(NULL, "-pc_hypre_pilut_factorrowsize", std::to_string(fill).c_str());
    //PCSetFromOptions(precon_);
}

int PetscSparseMatrixSolver::nIters() const
{
    PetscInt its;
    KSPGetIterationNumber(solver_, &its);
    return its;
}

Scalar PetscSparseMatrixSolver::error() const
{
    PetscReal err;
    KSPGetResidualNorm(solver_, &err);
    return err;
}

bool PetscSparseMatrixSolver::supportsMPI() const
{
    PCType type;
    PCGetType(precon_, &type);
    return strcmp(type, "hypre") || strcmp(type, "sor");
}

void PetscSparseMatrixSolver::printStatus(const std::string &msg) const
{
    comm_.printf("%s iterations = %d, error = %lf.\n", msg.c_str(), nIters(), error());
}
