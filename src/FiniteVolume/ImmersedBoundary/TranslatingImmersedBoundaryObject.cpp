#include "TranslatingImmersedBoundaryObject.h"


void TranslatingImmersedBoundaryObject::update(Scalar timeStep)
{
    shape() += timeStep*velocity_;
    updateCells();
}
