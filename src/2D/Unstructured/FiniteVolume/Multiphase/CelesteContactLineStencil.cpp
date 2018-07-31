#include "CelesteImmersedBoundary.h"
#include "Geometry/Intersection.h"

std::queue<Ref<const Cell>> CelesteImmersedBoundary::ContactLineStencil::cellQueue_;

std::unordered_set<Label> CelesteImmersedBoundary::ContactLineStencil::cellIdSet_;

CelesteImmersedBoundary::ContactLineStencil::ContactLineStencil(const Cell &cell,
                                                                const ImmersedBoundaryObject &ibObj,
                                                                Scalar theta, const ScalarFiniteVolumeField &gamma)
    :
      cell_(cell),
      link_(nullptr),
      theta_(theta)
{
    Vector2D ns = -ibObj.nearestEdgeUnitNormal(cell.centroid());

    //- The contact line normals
    //ncl1 = ns.rotate(-theta_);
    //ncl2 = ns.rotate(theta_);

    Ray2D r1 = Ray2D(cell_.centroid(), ns.rotate(M_PI_2 - theta_));
    Ray2D r2 = Ray2D(cell_.centroid(), ns.rotate(theta_ - M_PI_2));
    Point2D bp1 = ibObj.shape().intersections(r1)[0];
    Point2D bp2 = ibObj.shape().intersections(r2)[0];

    cellQueue_.push(std::cref(cell.grid().globalCells().nearestItem(bp1)));
    cellIdSet_.insert(cellQueue_.back().get().id());

    cellQueue_.push(std::cref(cell.grid().globalCells().nearestItem(bp2)));
    cellIdSet_.insert(cellQueue_.back().get().id());

    std::pair<Scalar, Scalar> minDistSqr;
    std::pair<const CellLink*, const CellLink*> links(nullptr, nullptr);
    std::pair<LineSegment2D, LineSegment2D> cls;

    while(!cellQueue_.empty())
    {
        const Cell &cell = cellQueue_.front();

        for(const CellLink &nb: cell.cellLinks())
        {
            if(cellIdSet_.find(nb.cell().id()) == cellIdSet_.end() && (!links.first || !links.second)) // don't add to queue if a result has been found
            {
                cellQueue_.push(std::cref(nb.cell()));
                cellIdSet_.insert(nb.cell().id());
            }

            if(!ibObj.isInIb(cell.centroid()) && !ibObj.isInIb(nb.cell().centroid()))
            {
                LineSegment2D ln(cell.centroid(), nb.cell().centroid());

                auto xc1 = intersection(r1, ln);
                auto xc2 = intersection(r2, ln);

                if(xc1.second)
                {
                    Scalar distSqr = (xc1.first - r1.x0()).magSqr();

                    if(!links.first || distSqr < minDistSqr.first)
                    {
                        minDistSqr.first = distSqr;
                        links.first = &nb;
                        cls.first = LineSegment2D(r1.x0(), xc1.first);
                    }
                }

                if(xc2.second)
                {
                    Scalar distSqr = (xc2.first - r2.x0()).magSqr();

                    if(!links.second || distSqr < minDistSqr.second)
                    {
                        minDistSqr.second = distSqr;
                        links.second = &nb;
                        cls.second = LineSegment2D(r2.x0(), xc2.first);
                    }
                }
            }
        }

        cellQueue_.pop();
    }

    cellIdSet_.clear();

    Scalar g1 = links.first->linearInterpolate(gamma, cls.first.ptB());
    Scalar g2 = links.second->linearInterpolate(gamma, cls.second.ptB());

    if(theta_ < M_PI_2 && g1 > g2 || theta_ > M_PI_2 && g1 < g2)
    {
        link_ = links.first;
        cl_ = cls.first;
        ncl_ = ns.rotate(-theta);
        gamma_ = g1;
    }
    else
    {
        link_ = links.second;
        cl_ = cls.second;
        ncl_ = ns.rotate(theta);
        gamma_ = g2;
    }
}
