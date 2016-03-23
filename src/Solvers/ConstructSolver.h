#ifndef CONSTRUCT_SOLVER_H
#define CONSTRUCT_SOLVER_H

#include <memory>

#include "Solver.h"
#include "FiniteVolumeGrid2D.h"
#include "Input.h"

std::shared_ptr<Solver> constructSolver(const FiniteVolumeGrid2D& grid, const Input& input);

#endif
