#ifndef EQUATION_H
#define EQUATION_H

#include "SparseMatrix.h"
#include "FiniteVolumeGrid2D.h"
#include "Term.h"
#include "DiffusionTerm.h"
#include "AdvectionTerm.h"
#include "TemporalTerm.h"

template<class T>
class Equation
{
public:

    Equation(const FiniteVolumeGrid2D& grid, T&field);

    Equation<T>& operator=(const Term& term);

    Scalar solve();

    Scalar error() const { return spMat_.error(); }
    int iterations() const { return spMat_.nIterations(); }

    const SparseMatrix& matrix() const { return spMat_; }

    void relax(Scalar relaxationFactor);

private:
    SparseMatrix spMat_;
    SparseVector x_, b_;
    const FiniteVolumeGrid2D& grid_;
    T& field_;
};

#include "Equation.tpp"

#endif
