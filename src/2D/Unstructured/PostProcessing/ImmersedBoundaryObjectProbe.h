#ifndef PHASE_IMMERSED_BOUNDARY_OBJECT_PROBE_H
#define PHASE_IMMERSED_BOUNDARY_OBJECT_PROBE_H

#include "PostProcessing.h"
#include "FiniteVolume/ImmersedBoundary/ImmersedBoundaryObject.h"

class ImmersedBoundaryObjectProbe : public PostProcessing::Object
{
public:

    ImmersedBoundaryObjectProbe(int fileWriteFreq,
                                const std::weak_ptr<const ImmersedBoundaryObject> &ibObj,
                                const std::weak_ptr<const ScalarFiniteVolumeField> &field,
                                const Vector2D& probePos);

    void compute(Scalar time, bool force = false) override;

protected:

    std::string getFilename() const;

    std::weak_ptr<const ImmersedBoundaryObject> ibObj_;

    std::weak_ptr<const ScalarFiniteVolumeField> field_;

    Point2D probePos_;
};

#endif
