#include "TranslatingImmersedBoundaryObject.h"


TranslatingImmersedBoundaryObject::TranslatingImmersedBoundaryObject(const std::string &name,
                                                                     const Vector2D &velocity,
                                                                     Scalar omega,
                                                                     Label id,
                                                                     FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryObject(name,
                               id,
                               grid),
        velocity_(velocity),
        omega_(omega)
{

}

void TranslatingImmersedBoundaryObject::update(Scalar timeStep)
{
    shape() += timeStep * velocity_;
    updateCells();
}
