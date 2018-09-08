#ifndef PHASE_CELESTE_IMMERSED_BOUNDARY_H
#define PHASE_CELESTE_IMMERSED_BOUNDARY_H

#include "Geometry/StaticPolyLine2D.h"

#include "Celeste.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundary.h"

class CelesteImmersedBoundary: public Celeste
{
public:

    CelesteImmersedBoundary(const Input &input,
                            const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                            const std::shared_ptr<CellGroup> &fluidCells,
                            const std::weak_ptr<const ImmersedBoundary> &ib);

    Scalar theta(const ImmersedBoundaryObject &ibObj) const;

    void computeContactLineExtension(ScalarFiniteVolumeField &gamma) const;

    FiniteVolumeEquation<Scalar> contactLineBcs(ScalarFiniteVolumeField &gamma, Scalar timeStep) const;

    void appyFluidForces(const ScalarFiniteVolumeField &rho,
                         const ScalarFiniteVolumeField &mu,
                         const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &p,
                         const ScalarFiniteVolumeField &gamma,
                         const Vector2D &g,
                         DirectForcingImmersedBoundary &ib) const;

protected:

    class ContactLineStencil
    {
    public:

        ContactLineStencil(const ImmersedBoundaryObject &ibObj,
                           const Point2D &pt,
                           Scalar theta,
                           const ScalarFiniteVolumeField &gamma);

        void init();

        const CellLink &link() const
        { return *link_; }

        const StaticPolyLine2D<3>& cl() const
        { return cl_; }

        const Vector2D& ncl() const
        { return ncl_; }

        Vector2D tcl() const
        { return (cl_[2] - cl_[0]).unitVec(); }

        Scalar gamma() const
        { return gamma_; }

        template<class T>
        T interpolate(const FiniteVolumeField<T> &field) const
        { return link_->linearInterpolate(field, cl_[2]); }

    protected:

        static std::queue<Ref<const Cell>> cellQueue_;

        static std::unordered_set<Label> cellIdSet_;

        void init(const Ray2D &r1, const Ray2D &r2, const ScalarFiniteVolumeField &gamma);

        std::pair<StaticPolyLine2D<3>, const CellLink*> findIntersectingCellLink(const Ray2D &r, const ImmersedBoundaryObject &ibObj);

        const ImmersedBoundaryObject &ibObj_;

        const CellLink* link_;

        Scalar theta_;

        Scalar gamma_;

        StaticPolyLine2D<3> cl_;

        Vector2D ncl_;
    };

    void computeInterfaceNormals() override;

    std::weak_ptr<const ImmersedBoundary> ib_;

    std::unordered_map<std::string, Scalar> ibContactAngles_;

    std::unordered_map<std::string, CellGroup> contactLineExensionCells_;
};

#endif
