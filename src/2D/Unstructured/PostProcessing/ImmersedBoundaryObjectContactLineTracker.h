#ifndef PHASE_IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H
#define PHASE_IMMERSED_BOUNDARY_OBJECT_CONTACT_LINE_TRACKER_H

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/ImmersedBoundary/ImmersedBoundary.h"

#include "PostProcessing.h"

class ImmersedBoundaryObjectContactLineTracker: public PostProcessing::Object
{
public:
    ImmersedBoundaryObjectContactLineTracker(int fileWriteFreq,
                                             const std::weak_ptr<const ScalarFiniteVolumeField> &gamma,
                                             const std::weak_ptr<const ImmersedBoundary> &ib);

    void compute(Scalar time);

private:

    std::weak_ptr<const ImmersedBoundary> ib_;

    std::weak_ptr<const ScalarFiniteVolumeField> gamma_;

};

#endif
