#ifndef GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H
#define GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H

#include <memory>

#include "PostProcessing.h"
#include "GhostCellImmersedBoundaryObject.h"

class GhostCellImmersedBoundaryObjectContactLineTracker: public PostProcessing
{
public:
    GhostCellImmersedBoundaryObjectContactLineTracker(const Solver& solver,
                                                      std::shared_ptr<ImmersedBoundaryObject> gIbObj);

    void compute(Scalar time);

private:
    std::shared_ptr<GhostCellImmersedBoundaryObject> gIbObj_;
};


#endif
