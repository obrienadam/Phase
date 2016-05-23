#ifndef SURFACE_TENSION_FORCE_H
#define SURFACE_TENSION_FORCE_H

#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "Input.h"

class SurfaceTensionForce
{
public:

    SurfaceTensionForce(const Input &input, const ScalarFiniteVolumeField& gamma);

    virtual VectorFiniteVolumeField compute() = 0;

protected:

    Scalar sigma_;
    const ScalarFiniteVolumeField &gamma_;
};

//- Header files for the available methods
#include "ContinuumSurfaceForce.h"
#include "SurfaceTensionForce.h"
#include "HeightFunction.h"

#endif
