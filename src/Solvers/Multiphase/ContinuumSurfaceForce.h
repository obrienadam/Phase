#ifndef CONTINUUM_SURFACE_FORCE_H
#define CONTINUUM_SURFACE_FORCE_H

#include "SurfaceTensionForce.h"
#include "CellSearch.h"

class ContinuumSurfaceForce : public SurfaceTensionForce
{
public:

    ContinuumSurfaceForce(const Input& input, const ScalarFiniteVolumeField &gamma);

    virtual VectorFiniteVolumeField compute();

private:

    void constructSmoothingKernels();
    void computeInterfaceNormals();
    void computeCurvature();
    Vector2D computeContactLineNormal(const Vector2D& gradGamma, const Vector2D& wallNormal);

    std::vector< std::vector< Ref<const Cell> > > cellRangeSearch_;
    Scalar kernelWidth_;

    ScalarFiniteVolumeField gammaTilde_, kappa_;
    VectorFiniteVolumeField n_;
};

#endif
