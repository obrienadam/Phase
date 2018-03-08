#include "DirectForcingImmersedBoundaryObject.h"

DirectForcingImmersedBoundaryObject::DirectForcingImmersedBoundaryObject(const std::string &name, Label id,
                                                                         const ImmersedBoundary &ib,
                                                                         const std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        ImmersedBoundaryObject(name, id, ib, grid)
{

}