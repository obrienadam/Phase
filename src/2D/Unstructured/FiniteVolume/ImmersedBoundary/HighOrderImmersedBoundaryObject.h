#ifndef PHASE_HIGH_ORDER_IMMERSED_BOUNDARY_OBJECT_H
#define PHASE_HIGH_ORDER_IMMERSED_BOUNDARY_OBJECT_H

#include "Math/Matrix.h"
#include "Math/StaticMatrix.h"

#include "ImmersedBoundaryObject.h"

class HighOrderImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    HighOrderImmersedBoundaryObject(const std::string &name,
                                    const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                    const std::shared_ptr<CellGroup> &solverCells);

    Type type() const
    { return HIGH_ORDER; }

    void updateCells();

    FiniteVolumeEquation<Scalar> bcs(ScalarFiniteVolumeField &phi) const;

    FiniteVolumeEquation<Vector2D> bcs(VectorFiniteVolumeField &u) const
    {}

    FiniteVolumeEquation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    FiniteVolumeEquation<Scalar> contactLineBcs(ScalarFiniteVolumeField &gamma, Scalar theta) const;

    void computeForce(Scalar rho,
                      Scalar mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const Vector2D &g = Vector2D(0., 0.));

private:

    //void constructDirichletCoeffsQuad();

    void constructDirichletCoeffs();

    //void constructNeumannCoeffs();

    void constructNeumannCoeffs();

    std::vector<std::vector<Ref<const Cell>>> stDCells_, stNCells_;
    std::vector<std::vector<Point2D>> bps_;
    std::vector<std::vector<Vector2D>> bns_;
    std::vector<StaticMatrix<6, 9>> Ad_, An_;
    std::vector<StaticMatrix<1, 9>> bd_, bn_;
};


#endif
