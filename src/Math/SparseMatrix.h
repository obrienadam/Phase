#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>

#include <eigen3/Eigen/Sparse>

#include "Types.h"
#include "SparseVector.h"

typedef Eigen::SparseMatrix<Scalar> SparseMatrix;
typedef Eigen::Triplet<Scalar> Triplet;

typedef Eigen::BiCGSTAB< SparseMatrix, Eigen::IncompleteLUT<Scalar> > BiCGSTABIncompleteLUT;
typedef Eigen::BiCGSTAB< SparseMatrix, Eigen::IdentityPreconditioner > BiCGSTABIdentityPreconditioner;
typedef Eigen::BiCGSTAB< SparseMatrix, Eigen::DiagonalPreconditioner<Scalar> > BiCGSTABDiagonalPreconditioner;

#endif
