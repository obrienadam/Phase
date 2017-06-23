#ifndef CONTINUUM_SURFACE_FORCE_H
#define CONTINUUM_SURFACE_FORCE_H

#include "SurfaceTensionForce.h"

class ContinuumSurfaceForce : public SurfaceTensionForce
{
public:

    ContinuumSurfaceForce(const Input &input,
                        const ScalarFiniteVolumeField& gamma,
                        const ScalarFiniteVolumeField& rho,
                        const ScalarFiniteVolumeField& mu,
                        const VectorFiniteVolumeField& u,
                        VectorFiniteVolumeField& gradGamma);

    virtual void compute();

    virtual const ScalarFiniteVolumeField& gammaTilde() const { return *gammaTilde_; }

    virtual const VectorFiniteVolumeField& gradGammaTilde() const { return *gradGammaTilde_; }

    virtual void registerSubFields(Solver& solver);

    void constructSmoothingKernels();

protected:

    virtual void computeGradGammaTilde();
    virtual void computeInterfaceNormals();
    virtual void computeCurvature();
    virtual void interpolateCurvatureFaces();

    std::vector< std::vector< Ref<const Cell> > > cellRangeSearch_;
    Scalar kernelWidth_;

    Scalar curvatureCutoffTolerance_;

    std::shared_ptr<ScalarFiniteVolumeField> gammaTilde_;
    std::shared_ptr<VectorFiniteVolumeField> gradGammaTilde_;
};



#endif
