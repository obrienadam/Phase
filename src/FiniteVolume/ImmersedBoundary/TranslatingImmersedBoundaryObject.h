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
    Vector2D velocity(const Point2D &point) const
    { return velocity_ + omega_ * (point - shapePtr_->centroid()).tangentVec(); }

    Vector2D trajectory() const
    { return velocity_; }

private:

    Scalar omega_ = 0.;
    Vector2D velocity_;

};

#endif
