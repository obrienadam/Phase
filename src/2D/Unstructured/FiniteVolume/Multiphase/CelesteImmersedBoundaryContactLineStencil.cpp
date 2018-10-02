#include "Geometry/Intersection.h"

#include "CelesteImmersedBoundary.h"

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

    if(theta_ == M_PI_2)
        init(Ray2D(pt, ns), gamma);
    else
        init(Ray2D(pt, ns.rotate(M_PI_2 - theta)), Ray2D(pt, theta_ - M_PI_2), gamma);
}

void CelesteImmersedBoundary::ContactLineStencil::init(const Ray2D &r, const ScalarFiniteVolumeField &gamma)
{
    auto l = findIntersectingCellLink(r, ibObj_);

    if(!l.second)
        return;

    cl_ = l.first;
    link_ = l.second;

    gamma_ = link_->linearInterpolate(gamma, l.first[2]);

    Vector2D nl = link_->rc().unitVec();
    nl = gamma(link_->self()) >= gamma(link_->cell()) ? nl : -nl;

    Vector2D ncl =r.r().rotate(M_PI_2);
    ncl_ = (dot(ncl, nl) * ncl).unitVec();
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

    //    Scalar dg = gamma(link_->cell()) - gamma(link_->self());
    //    if(dot(dg * link_->rCellVec(), ncl_) > 0.)
    //        ncl_ = -ncl_;
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
