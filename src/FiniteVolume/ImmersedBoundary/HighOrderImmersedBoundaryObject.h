#ifndef HIGH_ORDER_IMMERSED_BOUNDARY_OBJECT_H
#define HIGH_ORDER_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"
#include "Matrix.h"
#include "StaticMatrix.h"

class HighOrderImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    HighOrderImmersedBoundaryObject(const std::string &name,
                                    Label id,
                                    FiniteVolumeGrid2D &grid);

    Type type() const
    { return HIGH_ORDER; }

    void updateCells();

    Equation<Scalar> bcs(ScalarFiniteVolumeField &phi) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField &u) const
    {}

    Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    Equation<Scalar> contactLineBcs(ScalarFiniteVolumeField &gamma, Scalar theta) const;

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
