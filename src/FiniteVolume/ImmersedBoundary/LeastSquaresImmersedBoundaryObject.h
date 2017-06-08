#ifndef LEAST_SQUARES_IMMERSED_BOUNDARY_OBJECT
#define LEAST_SQUARES_IMMERSED_BOUNDARY_OBJECT

#include "ImmersedBoundaryObject.h"
#include "LeastSquaresStencil.h"

class LeastSquaresImmersedBoundaryObject: public ImmersedBoundaryObject
{
public:
    LeastSquaresImmersedBoundaryObject(const std::string &name,
                                       Label id,
                                       FiniteVolumeGrid2D &grid);

    void update(Scalar timeStep);

    void updateCells();

    Equation<Scalar> bcs(ScalarFiniteVolumeField& field) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField& field) const;

protected:

    void constructStencils();

    std::vector<LeastSquaresStencil> stencils_;
    bool updateCellsCalled_ = false;
};

#endif
