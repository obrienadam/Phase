#ifndef EQUATION_H
#define EQUATION_H

#include "SparseMatrix.h"
#include "FiniteVolumeGrid2D.h"
#include "Term.h"

template<class T>
class Equation
{
public:

    Equation(const FiniteVolumeGrid2D& grid, T&field);

    void operator==(const Term& term);

private:
    SparseMatrix spMat_;
    const FiniteVolumeGrid2D& grid_;
    T& field_;
};

#endif
