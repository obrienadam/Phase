#ifndef DIRECT_FORCING_IMMERSED_BOUNDARY_OBJECT_H
#define DIRECT_FORCING_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"

class DirectForcingImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:
    //- Constructors, one for circles, another for polygons
    DirectForcingImmersedBoundaryObject(const std::string &name,
                                        Label id,
                                        const ImmersedBoundary &ib,
                                        const std::shared_ptr<FiniteVolumeGrid2D> &grid);
};


#endif
