#include "OscillatingImmersedBoundaryObject.h"

OscillatingImmersedBoundaryObject::OscillatingImmersedBoundaryObject(const std::string &name,
                                                                     Vector2D freq,
                                                                     Vector2D amp,
                                                                     Vector2D phaseShift,
                                                                     Point2D origin,
                                                                     Label id,
                                                                     FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryObject(name, id, grid)
{
    amp_ = amp;
    freq_ = 2 * M_PI * freq;
    phaseShift_ = phaseShift * M_PI / 180.;
    origin_ = origin;
}

void OscillatingImmersedBoundaryObject::updateCells()
{
    ImmersedBoundaryObject::updateCells();
}

void OscillatingImmersedBoundaryObject::update(Scalar timeStep)
{
    shapePtr_->move(
            origin_ +
            Vector2D(amp_.x * sin(freq_.x * time_ + phaseShift_.x), amp_.y * sin(freq_.y * time_ + phaseShift_.y))
    );

    time_ += timeStep;
    updateCells();
}

Vector2D OscillatingImmersedBoundaryObject::acceleration() const
{
    return Vector2D(
            -amp_.x * pow(freq_.x, 2) * sin(freq_.x * time_ + phaseShift_.x),
            -amp_.y * pow(freq_.y, 2) * sin(freq_.y * time_ + phaseShift_.y)
    );
}

Vector2D OscillatingImmersedBoundaryObject::acceleration(const Point2D &point) const
{
    return acceleration();
}

Vector2D OscillatingImmersedBoundaryObject::velocity() const
{
    return Vector2D(
            amp_.x * freq_.x * cos(freq_.x * time_ + phaseShift_.x),
            amp_.y * freq_.y * cos(freq_.y * time_ + phaseShift_.y)
    );
}

Vector2D OscillatingImmersedBoundaryObject::velocity(const Point2D &point) const
{
    return velocity();
}
