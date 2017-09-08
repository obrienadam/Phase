#ifndef GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H
#define GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H

#include <memory>

#include "PostProcessingObject.h"
#include "GhostCellImmersedBoundaryObject.h"

class GhostCellImmersedBoundaryObjectContactLineTracker: public PostProcessingObject
{
public:
    GhostCellImmersedBoundaryObjectContactLineTracker(const Solver& solver);

    void compute(Scalar time);

private:
    std::vector<std::weak_ptr<GhostCellImmersedBoundaryObject>> gcIbObjs_;
};


#endif
