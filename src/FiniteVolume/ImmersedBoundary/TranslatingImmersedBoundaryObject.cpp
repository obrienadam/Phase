#include "TranslatingImmersedBoundaryObject.h"


TranslatingImmersedBoundaryObject::TranslatingImmersedBoundaryObject(const std::string &name, const Vector2D &velocity,
                                                                     Label id,
                                                                     FiniteVolumeGrid2D &grid)
    :
      ImmersedBoundaryObject(name,
                             id,
                             grid),
      velocity(velocity)
{

}

void TranslatingImmersedBoundaryObject::update(Scalar timeStep)
{
    shape() += timeStep*velocity;
    updateCells();
}
