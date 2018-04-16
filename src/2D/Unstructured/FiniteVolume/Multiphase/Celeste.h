#ifndef PHASE_CELESTE_H
#define PHASE_CELESTE_H

#include "System/Input.h"
#include "Math/Matrix.h"

#include "SurfaceTensionForce.h"

class Celeste : public SurfaceTensionForce
{
public:

    Celeste(const Input &input,
            const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
            const std::weak_ptr<ImmersedBoundary> &ib);

    void computeFaceInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma);

    void computeInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma);

protected:

    class CelesteStencil
    {
    public:

        CelesteStencil()
        {}

        CelesteStencil(const Cell& cell, bool weighted = false);

        CelesteStencil(const Cell& cell, const ImmersedBoundary& ib, bool weighted = false);

        void init(bool weighted = false);

        void init(const ImmersedBoundary& ib, bool weighted = false);

        void reset();

        const Cell& cell() const
        { return *cellPtr_; }

        bool weighted() const
        { return weighted_; }

        bool truncated() const
        { return truncated_; }

        Vector2D grad(const ScalarFiniteVolumeField& phi) const;

        Scalar div(const VectorFiniteVolumeField& u) const;

        Scalar kappa(const VectorFiniteVolumeField& n, const Celeste& fst) const;

        Scalar kappa(const VectorFiniteVolumeField& n, const ImmersedBoundary& ib, const Celeste& fst) const;

    private:

        void initMatrix();

        const Cell* cellPtr_ = nullptr;

        bool truncated_, weighted_;

        Matrix pInv_;

        std::vector<Ref<const Cell>> cells_;

        std::vector<Ref<const Face>> faces_;

        std::vector<std::pair<Ref<const Cell>, std::weak_ptr<const ImmersedBoundaryObject>>> compatPts_;
    };

    void computeGradGammaTilde(const ScalarFiniteVolumeField &gamma);

    void computeCurvature();

    void updateStencils();

    std::vector<CelesteStencil> kappaStencils_, gradGammaTildeStencils_;
};

#endif
