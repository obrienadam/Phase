#ifndef CONTINUUM_SURFACE_FORCE_H
#define CONTINUUM_SURFACE_FORCE_H

#include "SurfaceTensionForce.h"
#include "CellSearch.h"

class ContinuumSurfaceForce : public SurfaceTensionForce
{
public:

    ContinuumSurfaceForce(const Input& input, const ScalarFiniteVolumeField &gamma, const ScalarFiniteVolumeField& rho);

    virtual VectorFiniteVolumeField compute();
    const VectorFiniteVolumeField& n() const { return n_; }

    virtual const VectorFiniteVolumeField& gradGamma() const { return gradGammaTilde_; }

protected:

    void constructSmoothingKernels();
    void computeGradGammaTilde();
    virtual void computeInterfaceNormals();
    virtual void computeCurvature();

    std::vector< std::vector< Ref<const Cell> > > cellRangeSearch_;
    Scalar kernelWidth_;

    Scalar avgRho_;
    const ScalarFiniteVolumeField &rho_;
    VectorFiniteVolumeField gradGammaTilde_;
};

#endif
