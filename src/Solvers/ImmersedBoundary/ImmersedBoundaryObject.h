#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Circle.h"
#include "FiniteVolumeGrid2D.h"

class ImmersedBoundaryObject : public Circle
{
public:

    class ImmersedBoundaryStencil
    {
    public:

        ImmersedBoundaryStencil(const Cell& cell, const Point2D& bp, const FiniteVolumeGrid2D& grid);

        const Cell& cell() const { return cell_; }
        const Point2D& boundaryPoint() const { return bp_; }
        const Point2D& imagePoint() const { return ip_; }

    private:

        Ref<const Cell> cell_;
        Point2D bp_, ip_;
        std::vector< Ref<const Cell> > kNN_;
    };

    ImmersedBoundaryObject(const FiniteVolumeGrid2D& grid, const Point2D& center, Scalar radius);

    void constructStencils();
    const std::vector<ImmersedBoundaryStencil>& stencils() const { return ibStencils_; }

private:

    const FiniteVolumeGrid2D& grid_;
    std::vector<ImmersedBoundaryStencil> ibStencils_;
};

#endif
