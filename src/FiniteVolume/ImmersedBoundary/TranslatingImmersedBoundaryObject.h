#ifndef TRANSLATING_IMMERSED_BOUNDARY_OBJECT
#define TRANSLATING_IMMERSED_BOUNDARY_OBJECT

#include "ImmersedBoundaryObject.h"

class TranslatingImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    TranslatingImmersedBoundaryObject(const std::string &name,
                                      const Vector2D &velocity,
                                      Scalar omega,
                                      Label id,
                                      FiniteVolumeGrid2D &grid);

    void update(Scalar timeStep);

    //- Motion info
    Vector2D acceleration() const
    { return Vector2D(0., 0.); }

    Vector2D acceleration(const Point2D &point) const
    { return Vector2D(0., 0.); }

    Vector2D velocity() const
    { return velocity_; }

    Vector2D velocity(const Point2D &point) const
    { return velocity_; }

    Scalar angularVelocity() const
    { return 0.; }

private:

    Scalar omega_ = 0.;
    Vector2D velocity_;

};

#endif
