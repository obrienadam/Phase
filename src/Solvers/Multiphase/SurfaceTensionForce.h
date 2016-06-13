#ifndef SURFACE_TENSION_FORCE_H
#define SURFACE_TENSION_FORCE_H

#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "Input.h"

class SurfaceTensionForce
{
public:

    SurfaceTensionForce(const Input &input, const ScalarFiniteVolumeField& gamma, const VectorFiniteVolumeField& u);

    virtual VectorFiniteVolumeField compute() = 0;
    Vector2D computeContactLineNormal(const Vector2D& gradGamma, const Vector2D& wallNormal, const Vector2D &vel) const;

    const ScalarFiniteVolumeField& gamma() const { return gamma_; }
    virtual const ScalarFiniteVolumeField& gammaTilde() const { return gamma_; }
    virtual const VectorFiniteVolumeField& gradGamma() const = 0;

    const ScalarFiniteVolumeField& kappa() const { return kappa_; }
    const VectorFiniteVolumeField& n() const { return n_; }

    Scalar theta() const { return thetaAdv_; }
    Scalar sigma() const { return sigma_; }

protected:

    Scalar sigma_, thetaAdv_, thetaRec_;
    const ScalarFiniteVolumeField &gamma_;
    const VectorFiniteVolumeField &u_;
    ScalarFiniteVolumeField kappa_;
    VectorFiniteVolumeField n_;

    std::vector< Ref<const Patch> > contactAnglePatches_;
};

//- Header files for the available methods
#include "ContinuumSurfaceForce.h"
#include "SurfaceTensionForce.h"
#include "HeightFunction.h"

#endif
