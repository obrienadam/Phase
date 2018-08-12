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

    Ray2D r1 = Ray2D(cell_.centroid(), ns.rotate(M_PI_2 - theta_));
    Ray2D r2 = Ray2D(cell_.centroid(), ns.rotate(theta_ - M_PI_2));

    auto l1 = findIntersectingCellLink(r1, ibObj);
    auto l2 = findIntersectingCellLink(r2, ibObj);

    if(!l1.second)
    {
        throw Exception("CelesteImmersedBoundary::ContactLineStencil",
                        "ContactLineStencil",
                        "no intersection found. Cell id = "
                        + std::to_string(cell.centroid())
                        + "."
                        + " r = " + std::to_string(r1.x0()) + std::to_string(r1.r()));
    }
    else if(!l2.second)
    {
        throw Exception("CelesteImmersedBoundary::ContactLineStencil",
                        "ContactLineStencil",
                        "no intersection found. Cell id = "
                        + std::to_string(cell.centroid())
                        + "."
                        + " r = " + std::to_string(r1.x0()) + std::to_string(r1.r()));
    }

    Scalar g1 = l1.second->linearInterpolate(gamma, l1.first[2]);
    Scalar g2 = l2.second->linearInterpolate(gamma, l2.first[2]);

    if(theta_ < M_PI_2 && g1 > g2 || theta_ > M_PI_2 && g1 < g2)
    {
        link_ = l1.second;
        cl_ = l1.first;
        ncl_ = -ibObj.nearestEdgeUnitNormal(cl_[1]).rotate(-theta);
        gamma_ = g1;
    }
    else
    {
        link_ = l2.second;
        cl_ = l2.first;
        ncl_ = -ibObj.nearestEdgeUnitNormal(cl_[1]).rotate(theta);
        gamma_ = g2;
    }
}

std::pair<PolyLine2D, const CellLink*> CelesteImmersedBoundary::ContactLineStencil::findIntersectingCellLink(const Ray2D &r,
                                                                                                             const ImmersedBoundaryObject &ibObj)
{
    auto intersections = ibObj.shape().intersections(r);

    if(intersections.size() != 1)
        throw Exception("CelesteImmersedBoundary::ContactLineStencil",
                        "findIntersectingCellLink",
                        "no intersection found.");

    Vector2D bp = intersections.front();

    cellQueue_.push(std::cref(cell_.grid().globalCells().nearestItem(bp)));
    cellIdSet_.insert(cellQueue_.back().get().id());

    const CellLink *link = nullptr;
    Scalar minDistSqr;
    PolyLine2D cl;

    while(!cellQueue_.empty())
    {
        const Cell &cell = cellQueue_.front();
        cellQueue_.pop();

        for(const CellLink &nb: cell.cellLinks())
        {
            if(cellIdSet_.find(nb.cell().id()) == cellIdSet_.end() && !link) // don't add to queue if a result has been found
            {
                cellQueue_.push(std::cref(nb.cell()));
                cellIdSet_.insert(nb.cell().id());
            }

            if(ibObj.isInIb(cell) || ibObj.isInIb(nb.cell()))
                continue;

            LineSegment2D ln = LineSegment2D(cell.centroid(), nb.cell().centroid());
            auto xc = intersection(r, ln);

            if(xc.second)
            {
                Scalar distSqr = (xc.first - r.x0()).magSqr();

                if(!link || distSqr < minDistSqr)
                {
                    minDistSqr = distSqr;
                    link = &nb;
                    cl = {
                        r.x0(),
                        bp,
                        xc.first
                    };
                }
            }
        }
    }

    cellIdSet_.clear();

    return std::make_pair(cl, link);
}
