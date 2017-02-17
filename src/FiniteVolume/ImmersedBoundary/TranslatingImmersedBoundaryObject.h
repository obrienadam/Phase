#ifndef TRANSLATING_IMMERSED_BOUNDARY_OBJECT
#define TRANSLATING_IMMERSED_BOUNDARY_OBJECT

#include "ImmersedBoundaryObject.h"

class TranslatingImmersedBoundaryObject: public ImmersedBoundaryObject
{
public:

    TranslatingImmersedBoundaryObject(const std::string& name,
                                      const Vector2D& velocity,
                                      Label id,
                                      FiniteVolumeGrid2D& grid);

    void update(Scalar timeStep);

    Vector2D velocity;

private:

};

#endif
