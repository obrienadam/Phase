#include "CelesteImmersedBoundary.h"
#include "Geometry/Intersection.h"

std::queue<Ref<const Cell>> CelesteImmersedBoundary::ContactLineStencil::cellQueue_;

std::unordered_set<Label> CelesteImmersedBoundary::ContactLineStencil::cellIdSet_;

CelesteImmersedBoundary::ContactLineStencil::ContactLineStencil(const ImmersedBoundaryObject &ibObj,
                                                                const Point2D &pt,
                                                                Scalar theta,
                                                                const ScalarFiniteVolumeField &gamma)
    :
      ibObj_(ibObj),
      link_(nullptr),
      theta_(theta)
{
    Vector2D ns = -ibObj_.nearestEdgeUnitNormal(pt);
    Ray2D r1 = Ray2D(pt, ns.rotate(M_PI_2 - theta_));
    Ray2D r2 = Ray2D(pt, ns.rotate(theta_ - M_PI_2));

    init(r1, r2, gamma);
}

void CelesteImmersedBoundary::ContactLineStencil::init(const Ray2D &r1, const Ray2D &r2, const ScalarFiniteVolumeField &gamma)
{
    auto l1 = findIntersectingCellLink(r1, ibObj_);
    auto l2 = findIntersectingCellLink(r2, ibObj_);

    if(!l1.second || !l2.second)
    {
        throw Exception("CelesteImmersedBoundary::ContactLineStencil",
                        "ContactLineStencil",
                        "no gridline intersection found.");
    }

    Scalar g1 = l1.second->linearInterpolate(gamma, l1.first[2]);
    Scalar g2 = l2.second->linearInterpolate(gamma, l2.first[2]);

    if(theta_ < M_PI_2 && g1 > g2 || theta_ > M_PI_2 && g1 < g2)
    {
        link_ = l1.second;
        cl_ = l1.first;
        ncl_ = -ibObj_.nearestEdgeUnitNormal(cl_[1]).rotate(-theta_);
        gamma_ = g1;
    }
    else
    {
        link_ = l2.second;
        cl_ = l2.first;
        ncl_ = -ibObj_.nearestEdgeUnitNormal(cl_[1]).rotate(theta_);
        gamma_ = g2;
    }
}

std::pair<StaticPolyLine2D<3>, const CellLink *> CelesteImmersedBoundary::ContactLineStencil::findIntersectingCellLink(const Ray2D &r,
                                                                                                                       const ImmersedBoundaryObject &ibObj)
{
    auto intersections = ibObj.shape().intersections(r);

    if(intersections.size() < 1)
        intersections.push_back(ibObj.nearestIntersect(r.x0()));

//    if(intersections.size() != 1)
//        throw Exception("CelesteImmersedBoundary::ContactLineStencil",
//                        "findIntersectingCellLink",
//                        "contact line does not intersect boundary.");

    Vector2D bp = intersections.front();

    cellQueue_.push(std::cref(ibObj.ibCells().nearestItem(bp)));
    cellIdSet_.insert(cellQueue_.back().get().id());

    const CellLink *link = nullptr;
    Scalar minDistSqr;
    StaticPolyLine2D<3> cl;

    bool addToQueue = true;

    while(!cellQueue_.empty())
    {
        const Cell &cell = cellQueue_.front();
        cellQueue_.pop();

        for(const CellLink &nb: cell.cellLinks())
        {
            if(cellIdSet_.find(nb.cell().id()) == cellIdSet_.end() && addToQueue) // don't add to queue if a result has been found
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

        addToQueue = !link;
    }

    cellIdSet_.clear();

    return std::make_pair(cl, link);
}
