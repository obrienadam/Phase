#ifndef GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_FORCE_INTEGRATOR_H
#define GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_FORCE_INTEGRATOR_H

#include "PostProcessingObject.h"
#include "GhostCellImmersedBoundaryObject.h"

class GhostCellImmersedBoundaryObjectForceIntegrator: public PostProcessingObject
{
public:

    GhostCellImmersedBoundaryObjectForceIntegrator(const Solver& solver);

    void compute(Scalar time);

private:

    std::vector<std::weak_ptr<GhostCellImmersedBoundaryObject>> gcIbObjs_;

};


#endif
