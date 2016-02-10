#ifndef EQUATION_H
#define EQUATION_H

#include "SparseMatrix.h"
#include "FiniteVolumeGrid2D.h"

class Equation
{
public:

    Equation(const FiniteVolumeGrid2D& grid);

private:
    SparseMatrix spMat_;
    const FiniteVolumeGrid2D& grid_;
};

#endif
