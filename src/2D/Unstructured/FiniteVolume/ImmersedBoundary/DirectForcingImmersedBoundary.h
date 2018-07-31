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

    void computeForcingTerm(const VectorFiniteVolumeField &u, Scalar timeStep, VectorFiniteVolumeField &fib) const;

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

    virtual void computeForce(Scalar rho,
                              Scalar mu,
                              const VectorFiniteVolumeField &u,
                              const ScalarFiniteVolumeField &p,
                              const Vector2D &g = Vector2D(0., 0.));

//    virtual void computeForce(const ScalarFiniteVolumeField &rho,
//                              const ScalarFiniteVolumeField &mu,
//                              const VectorFiniteVolumeField &u,
//                              const ScalarFiniteVolumeField &p,
//                              const ScalarFiniteVolumeField &gamma,
//                              const SurfaceTensionForce &ft,
//                              const Vector2D &g = Vector2D(0., 0.));


private:

    CellGroup ibCells_, solidCells_;

};


#endif
