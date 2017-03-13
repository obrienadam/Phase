#include "OscillatingImmersedBoundaryObject.h"

OscillatingImmersedBoundaryObject::OscillatingImmersedBoundaryObject(const std::string &name,
                                                                     Scalar amplitude,
                                                                     Scalar frequency,
                                                                     Label id,
                                                                     FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryObject(name, id, grid)
{
    amp_ = amplitude;
    omega_ = 2. * M_PI * frequency;
}

void OscillatingImmersedBoundaryObject::update(Scalar timeStep)
{
    if (time_ == 0.)
        origin_ = shape().centroid();

    Vector2D pos = origin_ + amp_ * sin(omega_ * (time_ += timeStep));
    shape() += pos - shape().centroid();

    updateCells();
}
