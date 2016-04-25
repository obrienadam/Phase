#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include "Multiphase.h"
#include "Circle.h"

class ImmersedBoundary : public Multiphase
{
public:

    ImmersedBoundary(const FiniteVolumeGrid2D& grid, const Input& input);

protected:
    std::vector<Circle> ibObjects_;
    std::vector< std::vector< Ref<const Cell> > > kNN_;
};

#endif
