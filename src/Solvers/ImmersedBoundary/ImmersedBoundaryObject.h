#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Circle.h"
#include "FiniteVolumeGrid2D.h"
#include "CellSearch.h"
#include "ContinuumSurfaceForce.h"

class ImmersedBoundaryObject : public Circle
{
public:

    enum BoundaryType{FIXED, NORMAL_GRADIENT, CONTACT_ANGLE};

    class ImmersedBoundaryStencil
    {
    public:

        ImmersedBoundaryStencil(const Cell& cell,
                                const Point2D& bp,
                                const FiniteVolumeGrid2D& grid,
                                const CellSearch &cs);

        const Cell& cell() const { return cell_; }
        const std::vector< Ref<const Cell> >& kNearestNeighbours() const { return kNN_; }

        const Point2D& boundaryPoint() const { return bp_; }
        const Point2D& imagePoint() const { return ip_; }

    private:

        Ref<const Cell> cell_;
        std::vector< Ref<const Cell> > kNN_;
        Point2D bp_, ip_;
    };

    ImmersedBoundaryObject(const FiniteVolumeGrid2D& grid, const std::shared_ptr<SurfaceTensionForce> &csfPtr);
    ImmersedBoundaryObject(const FiniteVolumeGrid2D &grid, const std::shared_ptr<SurfaceTensionForce> &csfPtr, const Point2D &center, Scalar radius);

    void init(const Point2D& center, Scalar radius);

    const FiniteVolumeGrid2D& grid() const { return grid_; }

    void constructStencils();

    const std::vector<ImmersedBoundaryStencil>& stencils() const { return ibStencils_; }

    void addBoundaryType(const std::string& fieldName, BoundaryType boundaryType);
    BoundaryType boundaryType(const std::string& fieldName) const { return boundaryTypes_.find(fieldName)->second; }

    const SurfaceTensionForce& csf() const { return *csf_; }

private:

    const FiniteVolumeGrid2D &grid_;
    std::shared_ptr<SurfaceTensionForce> csf_;
    std::vector<ImmersedBoundaryStencil> ibStencils_;

    std::map<std::string, BoundaryType> boundaryTypes_;
};

#endif
