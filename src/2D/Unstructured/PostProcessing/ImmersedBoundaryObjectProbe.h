#ifndef PHASE_IMMERSED_BOUNDARY_OBJECT_PROBE_H
#define PHASE_IMMERSED_BOUNDARY_OBJECT_PROBE_H

#include "PostProcessing.h"

class ImmersedBoundaryObjectProbe : public PostProcessing::Object
{
public:

    ImmersedBoundaryObjectProbe(const Solver &solver,
                                const std::string& ibObjName,
                                const std::string& fieldName,
                                const Vector2D& probePos);

    void compute(Scalar time);

protected:

    std::string getFilename() const;

    std::weak_ptr<const ImmersedBoundaryObject> ibObj_;

    std::weak_ptr<const ScalarFiniteVolumeField> field_;

    Point2D probePos_;
};

#endif
