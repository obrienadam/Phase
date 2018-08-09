#ifndef PHASE_BOUNDARY_PATCH_H
#define PHASE_BOUNDARY_PATCH_H

#include "Set.h"
#include "Face.h"

class BoundaryPatch: public Set<Face>
{
public:

    BoundaryPatch(const std::string &name);

    bool add(const Face& face) override;

protected:

};

#endif
