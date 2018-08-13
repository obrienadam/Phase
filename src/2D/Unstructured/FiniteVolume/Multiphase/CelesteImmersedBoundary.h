#ifndef PHASE_CELESTE_IMMERSED_BOUNDARY_H
#define PHASE_CELESTE_IMMERSED_BOUNDARY_H

#include "Geometry/PolyLine2D.h"

#include "Celeste.h"
#include "FiniteVolume/ImmersedBoundary/ImmersedBoundary.h"

class CelesteImmersedBoundary: public Celeste
{
public:

    CelesteImmersedBoundary(const Input &input,
                            const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                            const std::shared_ptr<CellGroup> &fluidCells,
                            const std::weak_ptr<const ImmersedBoundary> &ib);

    Scalar theta(const ImmersedBoundaryObject &ibObj) const;

    void computeContactLineExtension(ScalarFiniteVolumeField &gamma) const;

    void contactLineBcs(FiniteVolumeEquation<Scalar> &gammaEqn) const;

protected:

    class ContactLineStencil
    {
    public:

        ContactLineStencil(const Cell &cell,
                           const ImmersedBoundaryObject &ibObj,
                           Scalar theta,
                           const ScalarFiniteVolumeField &gamma);

        void init();

        const CellLink &link() const
        { return *link_; }

        const PolyLine2D& cl() const
        { return cl_; }

        const Vector2D& ncl() const
        { return ncl_; }

        Scalar gamma() const
        { return gamma_; }

    protected:

        std::pair<PolyLine2D, const CellLink*> findIntersectingCellLink(const Ray2D &r, const ImmersedBoundaryObject &ibObj);

        static std::queue<Ref<const Cell>> cellQueue_;

        static std::unordered_set<Label> cellIdSet_;

        const Cell& cell_;

        const CellLink* link_;


        Scalar theta_;

        Scalar gamma_;

        PolyLine2D cl_;

        Vector2D ncl_;
    };

    void computeInterfaceNormals() override;

    void computeCurvature() override;

    std::weak_ptr<const ImmersedBoundary> ib_;

    std::unordered_map<std::string, Scalar> ibContactAngles_;

    std::unordered_map<std::string, CellGroup> contactLineExensionCells_;
};

#endif
