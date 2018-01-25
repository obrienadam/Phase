#ifndef IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H
#define IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H

#include "PostProcessingObject.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class ImmersedBoundaryObjectContactLineTracker: public PostProcessingObject
{
public:
    ImmersedBoundaryObjectContactLineTracker(const Solver& solver);

    void compute(Scalar time);

private:


};

#endif