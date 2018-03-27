#ifndef PHASE_IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H
#define PHASE_IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"

#include "PostProcessing.h"

class ImmersedBoundaryObjectContactLineTracker: public PostProcessing::Object
{
public:
    ImmersedBoundaryObjectContactLineTracker(const Solver& solver);

    void compute(Scalar time);

private:

    std::weak_ptr<const ImmersedBoundary> ib_;

    std::weak_ptr<const ScalarFiniteVolumeField> gamma_;

};

#endif