#ifndef EULER_LAGRANGE_IMMERSEDBOUNDARY_OBJECT_H
#define EULER_LAGRANGE_IMMERSEDBOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"
#include "NotImplementedException.h"

class EulerLagrangeImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:
    EulerLagrangeImmersedBoundaryObject(const std::string &name,
                                        Label id,
                                        const ImmersedBoundary &ib,
                                        const std::shared_ptr<FiniteVolumeGrid2D> &grid);

    Type type() const
    { return EULER_LAGRANGE; }

    void initCircle(const Point2D &center, Scalar radius);

    void updateCells();

    //- These methods aren't used, this is strictly a corrective method
    Equation<Scalar> bcs(ScalarFiniteVolumeField &field) const
    { throw NotImplementedException("EulerLagrangeImmersedBoundaryObject", "bcs"); }

    Equation<Vector2D> bcs(VectorFiniteVolumeField &field) const
    { throw NotImplementedException("EulerLagrangeImmersedBoundaryObject", "bcs"); }

    Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    void computeForce(Scalar rho,
                              Scalar mu,
                              const VectorFiniteVolumeField &u,
                              const ScalarFiniteVolumeField &p,
                              const Vector2D &g = Vector2D(0., 0.))
    {}

    void correctVelocity(VectorFiniteVolumeField &u) const;

    Scalar kernel(const Point2D &x, const Point2D &xl) const;

private:

    void initLagrangePoints(int nLagrangePoints);

    void updateLagrangePoints()
    { initLagrangePoints(lagrangePoints_.size()); }

    std::shared_ptr<TrilinosAmesosSparseMatrixSolver> solver_;

    Scalar h_;

    std::vector<Point2D> lagrangePoints_;

    std::vector<std::vector<Ref<const Cell>>> lagrangeStencils_;

};


#endif
