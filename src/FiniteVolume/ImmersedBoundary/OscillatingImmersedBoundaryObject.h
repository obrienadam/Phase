#ifndef OSCILLATING_IMMERSED_BOUNDARY_OBJECT
#define OSCILLATING_IMMERSED_BOUNDARY_OBJECT

#include "ImmersedBoundaryObject.h"

class OscillatingImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    OscillatingImmersedBoundaryObject(const std::string &name,
                                      Vector2D amp,
                                      Vector2D freq,
                                      Vector2D phaseShift,
                                      Point2D origin,
                                      Label id,
                                      FiniteVolumeGrid2D &grid);

    virtual void updateCells();

    void update(Scalar timeStep);

    //- Motion info

    Vector2D acceleration() const;

    Vector2D acceleration(const Point2D &point) const;

    Vector2D velocity() const;

    Vector2D velocity(const Point2D &point) const;

private:

    Scalar time_ = 0.;
    Vector2D amp_, freq_, phaseShift_;
    Vector2D origin_;
    Vector2D prevVelocity_;

};

#endif
