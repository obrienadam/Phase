#ifndef PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_H
#define PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_H

#include "ImmersedBoundary.h"
#include "Math/StaticMatrix.h"
#include "Math/Matrix.h"

class DirectForcingImmersedBoundary : public ImmersedBoundary
{
public:

    class LeastSquaresQuadraticStencil;

    //- Constructors, one for circles, another for polygons
    DirectForcingImmersedBoundary(const Input &input,
                                  const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                  const std::shared_ptr<CellGroup> &domainCells);

    void updateCells();

    FiniteVolumeEquation<Vector2D> computeForcingTerm(const VectorFiniteVolumeField &u,
                                                      Scalar timeStep,
                                                      VectorFiniteVolumeField &fib) const;

    void computeFaceForcingTerm(const VectorFiniteVolumeField &u,
                                Scalar timeStep,
                                VectorFiniteVolumeField &fib) const;

    void computeForcingTerm(const ScalarFiniteVolumeField &rho,
                            const VectorFiniteVolumeField &u,
                            Scalar timeStep,
                            VectorFiniteVolumeField &fib) const;

    virtual FiniteVolumeEquation<Scalar> bcs(ScalarFiniteVolumeField &phi) const
    { return FiniteVolumeEquation<Scalar>(phi); }

    virtual FiniteVolumeEquation<Vector2D> bcs(VectorFiniteVolumeField &phi) const
    { return FiniteVolumeEquation<Vector2D>(phi); }

    virtual FiniteVolumeEquation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const
    { return FiniteVolumeEquation<Vector2D>(u); }

    virtual void applyHydrodynamicForce(Scalar rho,
                                        Scalar mu,
                                        const VectorFiniteVolumeField &u,
                                        const ScalarFiniteVolumeField &p,
                                        const Vector2D &g = Vector2D(0., 0.)) override;

    virtual void applyHydrodynamicForce(const ScalarFiniteVolumeField &rho,
                                        const ScalarFiniteVolumeField &mu,
                                        const VectorFiniteVolumeField &u,
                                        const ScalarFiniteVolumeField &p,
                                        const Vector2D &g = Vector2D(0., 0.)) override;

    void applyHydrodynamicForce(const VectorFiniteVolumeField &fib);

    //- Cell group access

    const CellGroup &localIbCells() const
    { return localIbCells_; }

    const CellGroup &globalIbCells() const
    { return globalIbCells_; }

    const CellGroup &localSolidCells() const
    { return localSolidCells_; }

    const CellGroup &globalSolidCells() const
    { return globalSolidCells_; }

private:

    CellGroup localIbCells_, globalIbCells_;

    CellGroup localSolidCells_, globalSolidCells_;
};


#endif
