#include "Geometry/Intersection.h"

#include "CelesteImmersedBoundary.h"

std::queue<Ref<const Cell>> CelesteImmersedBoundary::ContactLineStencil::cellQueue_;

std::unordered_set<Label> CelesteImmersedBoundary::ContactLineStencil::cellIdSet_;

CelesteImmersedBoundary::ContactLineStencil::ContactLineStencil(const ImmersedBoundaryObject &ibObj,
                                                                const Point2D &pt,
                                                                Scalar theta,
                                                                const ScalarFiniteVolumeField &gamma)
    :
      ibObj_(&ibObj),
      cellA_(nullptr),
      cellB_(nullptr),
      face_(nullptr),
      theta_(theta)
{
    cl_[0] = pt;
    ns_ = -ibObj_->nearestEdgeUnitNormal(cl_[0]);

    if(theta_ == M_PI_2)
        init(gamma);
    else
        init(Ray2D(cl_[0], ns_.rotate(M_PI_2 - theta_)), Ray2D(cl_[0], ns_.rotate(theta_ - M_PI_2)), gamma);
}

void CelesteImmersedBoundary::ContactLineStencil::init(const ScalarFiniteVolumeField &gamma)
{
    findStencilCells(Ray2D(cl_[0], ns_));

    if(!isValid())
    {
        std::ostringstream os;

        os << "Failed to find any suitable stencil." << "\n"
           << "Start point = " << cl_[0] << ", ns = "  << ns_ << ".\n"
           << "Cells searched = " << cellIdSet_.size() << ".";

        throw Exception("CelesteImmersedBoundary::ContactLineStencil",
                        "init",
                        os.str());
    }

    gamma_ = interpolate(gamma);

    Vector2D nl;

    if(cellA_ && cellB_)
    {
        nl = (cellB_->centroid() - cellA_->centroid()).unitVec();
        nl = gamma(*cellA_) >= gamma(*cellB_) ? nl : -nl;
    }
    else if(cellA_ && face_)
    {
        nl = (face_->centroid() - cellA_->centroid()).unitVec();
        nl = gamma(*cellA_) >= gamma(*face_) ? nl : -nl;
    }

    ncl_ = ns_.rotate(M_PI_2);
    ncl_ = (dot(ncl_, nl) * ncl_).unitVec();

    if(std::isnan(ncl_.x) || std::isnan(ncl_.y))
        ncl_ = ns_.rotate(M_PI_2);
}

void CelesteImmersedBoundary::ContactLineStencil::init(const Ray2D &r1, const Ray2D &r2, const ScalarFiniteVolumeField &gamma)
{
    findStencilCells(r1);
    ContactLineStencil c1 = *this;

    if(!c1.isValid())
    {
        findStencilCells(Ray2D(cl_[0], ns_));
        c1 = *this;
    }

    findStencilCells(r2);
    ContactLineStencil c2 = *this;

    if(!c2.isValid())
    {
        findStencilCells(Ray2D(cl_[0], ns_));
        c2 = *this;
    }

    if(c1.isValid() && c2.isValid())
    {
        Scalar g1 = c1.interpolate(gamma);
        Scalar g2 = c2.interpolate(gamma);

        if(theta_ < M_PI_2 && g1 > g2 || theta_ > M_PI_2 && g1 < g2)
        {
            *this = c1;
            gamma_ = g1;
            ncl_ = ns_.rotate(-theta_);
        }
        else
        {
            *this = c2;
            gamma_ = g2;
            ncl_ = ns_.rotate(theta_);
        }
    }
    else
    {
        std::cout << "Warning, failed to find valid contact line stencil on proc " << gamma.grid()->comm().rank()
                  << ". Initiating normal stencil stencil...\n";
        init(gamma);
    }
}

void CelesteImmersedBoundary::ContactLineStencil::findStencilCells(const Ray2D &r)
{
    auto intersections = ibObj_->shape().intersections(r);

    if(intersections.empty())
        intersections.emplace_back(r.x0());

    if(intersections.empty())
        throw Exception("CelesteImmersedBoundary::ContactLineStencil",
                        "findStencilCells",
                        "contact line does not intersect boundary.");

    Scalar minDistSqr;
    bool foundIntersection = false;

    cellA_ = nullptr;
    cellB_ = nullptr;
    face_ = nullptr;

    Vector2D bp = intersections.front();
    cellQueue_.emplace(ibObj_->cells().nearestItem(bp));
    cellIdSet_.clear();
    cellIdSet_.emplace(cellQueue_.back().get().id());

    int searches = 0;

    while(!cellQueue_.empty())
    {
        const Cell &cell = cellQueue_.front();
        cellQueue_.pop();

        for(const CellLink &nb: cell.cellLinks())
        {
            if(cellIdSet_.find(nb.cell().id()) == cellIdSet_.end() && !foundIntersection) // don't add to queue if a result has been found
            {
                cellQueue_.emplace(nb.cell());
                cellIdSet_.emplace(nb.cell().id());
            }

            if(ibObj_->isInIb(cell) || ibObj_->isInIb(nb.cell()))
                continue;

            LineSegment2D ln = LineSegment2D(cell.centroid(), nb.cell().centroid());
            auto xc = intersection(r, ln);

            if(xc.second)
            {
                Scalar distSqr = (xc.first - r.x0()).magSqr();

                if(!foundIntersection || distSqr < minDistSqr)
                {
                    foundIntersection = true;
                    minDistSqr = distSqr;

                    cellA_ = &nb.self();
                    cellB_ = &nb.cell();
                    face_ = nullptr;

                    cl_ = {
                        r.x0(),
                        bp,
                        xc.first
                    };
                }
            }
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            if(ibObj_->isInIb(cell) || ibObj_->isInIb(bd.face().centroid()))
                continue;

            LineSegment2D ln = LineSegment2D(cell.centroid(), bd.face().centroid());
            auto xc = intersection(r, ln);

            if(xc.second)
            {
                Scalar distSqr = (xc.first - r.x0()).magSqr();

                if(!foundIntersection || distSqr < minDistSqr)
                {
                    foundIntersection = true;
                    minDistSqr = distSqr;

                    cellA_ = &bd.self();
                    cellB_ = nullptr;
                    face_ = &bd.face();

                    cl_ = {
                        r.x0(),
                        bp,
                        xc.first
                    };
                }
            }
        }

        if(++searches > 144)
            break;
    }

    if(foundIntersection)
    {
        if(cellA_ && cellB_)
        {
            Scalar l1 = (cl_[2] - cellA_->centroid()).mag();
            Scalar l2 = (cl_[2] - cellB_->centroid()).mag();
            alpha_ = l2 / (l1 + l2);
        }
        else if(cellA_ && face_)
        {
            Scalar l1 = (cl_[2] - cellA_->centroid()).mag();
            Scalar l2 = (cl_[2] - face_->centroid()).mag();
            alpha_ = l2 / (l1 + l2);
        }
    }
}
