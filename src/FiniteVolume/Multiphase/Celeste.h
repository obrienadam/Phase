#ifndef CELESTE_H
#define CELESTE_H

#include "SurfaceTensionForce.h"
#include "Matrix.h"
#include "GhostCellStencil.h"

class Celeste : public SurfaceTensionForce
{
public:

    Celeste(const Input &input,
            const std::weak_ptr<ImmersedBoundary> &ib,
            ScalarFiniteVolumeField &gamma,
            const ScalarGradient &gradGamma,
            const ScalarFiniteVolumeField &rho,
            const ScalarFiniteVolumeField &mu,
            const VectorFiniteVolumeField &u);

    void computeFaces();

    void compute();

    void constructMatrices();

protected:

    class CelesteStencil
    {
    public:

        CelesteStencil() {}

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

    virtual void computeGradGammaTilde();

    virtual void computeCurvature();

    void computeCurvature(const ImmersedBoundary &ib);

    void updateStencils(const ImmersedBoundary& ib);

    std::vector<bool> modifiedStencil_;
    std::vector<CelesteStencil> kappaStencils_, gradGammaTildeStencils_;
};

#endif
