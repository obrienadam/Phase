#ifndef CONTINUUM_SURFACE_FORCE_H
#define CONTINUUM_SURFACE_FORCE_H

#include "SurfaceTensionForce.h"

class ContinuumSurfaceForce : public SurfaceTensionForce
{
public:

    ContinuumSurfaceForce(const Input& input,
                          const ScalarFiniteVolumeField& gamma,
                          const VectorFiniteVolumeField &u,
                          std::map<std::string, ScalarFiniteVolumeField>& fields);

    virtual VectorFiniteVolumeField compute();
    const VectorFiniteVolumeField& n() const { return n_; }

    virtual const ScalarFiniteVolumeField& gammaTilde() const { return gammaTilde_; }
    virtual const VectorFiniteVolumeField& gradGamma() const { return gradGammaTilde_; }

    void constructSmoothingKernels();

protected:

    void computeGradGammaTilde();
    virtual void computeInterfaceNormals();
    virtual void computeCurvature();

    std::vector< std::vector< Ref<const Cell> > > cellRangeSearch_;
    Scalar kernelWidth_;

    Scalar avgRho_;

    ScalarFiniteVolumeField &gammaTilde_;
    VectorFiniteVolumeField gradGammaTilde_;
};

#endif
