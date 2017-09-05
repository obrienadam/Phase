#ifndef CONSTRUCT_GRID_H
#define CONSTRUCT_GRID_H

#include <memory>

#include "Input.h"
#include "FiniteVolumeGrid2D.h"

std::shared_ptr<FiniteVolumeGrid2D> constructGrid(const Input &input, std::shared_ptr<Communicator> comm);

#endif
