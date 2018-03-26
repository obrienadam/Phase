//
// Created by aobrien on 26/03/18.
//

#ifndef IMMERSED_BOUNDARY_OBJECT_PROBE_H
#define IMMERSED_BOUNDARY_OBJECT_PROBE_H

#include "PostProcessingObject.h"

class ImmersedBoundaryObjectProbe : public PostProcessingObject
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
