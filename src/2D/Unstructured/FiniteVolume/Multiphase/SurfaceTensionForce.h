#ifndef SURFACE_TENSION_FORCE_H
#define SURFACE_TENSION_FORCE_H

#include "FiniteVolume/Field/VectorFiniteVolumeField.h"
#include "FiniteVolume/Field/ScalarGradient.h"
#include "FiniteVolume/Equation/Equation.h"

class ImmersedBoundary;
class ImmersedBoundaryObject;

class SurfaceTensionForce : public VectorFiniteVolumeField
{
public:

    //- Constructor
    SurfaceTensionForce(const Input &input,
                        const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                        const std::weak_ptr<ImmersedBoundary> &ib);

    //- Internal field pointers
    const std::shared_ptr<ScalarFiniteVolumeField> &kappa() const
    { return kappa_; }

    const std::shared_ptr<ScalarFiniteVolumeField> &gammaTilde() const
    { return gammaTilde_; }

    const std::shared_ptr<VectorFiniteVolumeField> &n() const
    { return n_; }

    const std::shared_ptr<ScalarGradient> &gradGammaTilde() const
    { return gradGammaTilde_; }

    //- Return properties
    Scalar sigma() const
    { return sigma_; }

    Scalar theta(const FaceGroup &patch) const;

    Scalar theta(const ImmersedBoundaryObject &ibObj) const;

    //- Compute
    virtual void computeFaceInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma) = 0;

    virtual void computeInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma) = 0;

    virtual void computeInterfaceNormals();

    Vector2D contactLineNormal(const Vector2D &n,
                               const Point2D &pt,
                               const ImmersedBoundaryObject &ibObj) const;

    Vector2D contactLineNormal(const Cell &cell,
                               const Point2D &pt,
                               const ImmersedBoundaryObject &ibObj) const;

    Vector2D contactLineTangent(const Vector2D &n,
                               const Point2D &pt,
                               const ImmersedBoundaryObject &ibObj) const;

    Vector2D contactLineTangent(const Cell &cell,
                                const Point2D &pt,
                                const ImmersedBoundaryObject &ibObj) const;

    void smoothGammaField(const ScalarFiniteVolumeField &gamma);

    virtual Vector2D computeCapillaryForce(const ScalarFiniteVolumeField &gamma,
                                           const ImmersedBoundaryObject &ibObj) const;

    //- Misc special gamma boundary equations
    virtual Equation<Scalar> contactLineBcs(ScalarFiniteVolumeField &gamma);

protected:

    Scalar sigma_, kernelWidth_, eps_ = 1e-8;

    std::unordered_map<std::string, Scalar> ibContactAngles_, patchContactAngles_;

    std::weak_ptr<ImmersedBoundary> ib_;

    //- Fields, can share ownership
    std::shared_ptr<ScalarFiniteVolumeField> kappa_, gammaTilde_;

    std::shared_ptr<VectorFiniteVolumeField> n_;

    std::shared_ptr<ScalarGradient> gradGammaTilde_;
};

#endif
