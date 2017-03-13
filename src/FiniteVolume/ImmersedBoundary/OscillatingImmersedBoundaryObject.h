#ifndef OSCILLATING_IMMERSED_BOUNDARY_OBJECT
#define OSCILLATING_IMMERSED_BOUNDARY_OBJECT

#include "ImmersedBoundaryObject.h"

class OscillatingImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    OscillatingImmersedBoundaryObject(const std::string &name,
                                      Scalar amplitude,
                                      Scalar frequency,
                                      Label id,
                                      FiniteVolumeGrid2D &grid);

    void update(Scalar timeStep);

    //- Motion info
    Vector2D velocity(const Point2D &point) const
    { return amp_ * omega_ * cos(omega_ * time_); }

    Vector2D trajectory() const
    { return amp_ * omega_ * cos(omega_ * time_); }

private:

    Scalar time_ = 0.;
    Scalar amp_, omega_;
    Vector2D origin_;

};

#endif
