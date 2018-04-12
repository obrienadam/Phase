#ifndef PHASE_GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H
#define PHASE_GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H

#include "FiniteVolumeGrid2D/BilinearInterpolator.h"

#include "ImmersedBoundaryObject.h"

class GhostCellImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    class Stencil
    {
    public:

        Stencil(const Cell &cell,
                const ImmersedBoundaryObject &ibObj);

        Stencil(const Cell &cell,
                const ImmersedBoundaryObject &ibObj,
                const Vector2D &r);

        virtual Scalar bpValue(const ScalarFiniteVolumeField &phi) const = 0;

        virtual void init() = 0;

        const Cell &cell() const
        { return cell_; }

        const std::vector<Ref<const Cell>> &cells() const
        { return cells_; }

        const std::vector<Scalar> &coeffs() const
        { return coeffs_; }

        const Point2D &ip() const
        { return ip_; }

        const Point2D &bp() const
        { return bp_; }

        const Vector2D &nw() const
        { return nw_; }

        Scalar length() const
        { return (ip_ - cell_.get().centroid()).mag(); }

    protected:

        Ref<const Cell> cell_;

        std::vector<Ref<const Cell>> cells_;

        std::vector<Scalar> coeffs_;

        Point2D ip_, bp_, nw_;
    };

    class FixedStencil : public Stencil
    {
    public:

        FixedStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj)
                : Stencil(cell, ibObj), bi_(ibObj.grid())
        { init(); }

        void init();

        Scalar bpValue(const ScalarFiniteVolumeField &phi) const
        { return (bi_(phi) + phi(cell_)) / 2.; }

        Tensor2D bpGrad(const VectorFiniteVolumeField &u) const
        { return bi_.grad(u); }

    protected:

        BilinearInterpolator bi_;
    };

    class ZeroGradientStencil : public Stencil
    {
    public:

        ZeroGradientStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj)
                :
                Stencil(cell, ibObj),
                bi_(ibObj.grid())
        {
            d_ = -nw_.unitVec();
            init();
        }

        ZeroGradientStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj, const Vector2D &r)
                :
                Stencil(cell, ibObj, r),
                bi_(ibObj.grid())
        {
            d_ = -r.unitVec();
            init();
        }

        Scalar bpValue(const ScalarFiniteVolumeField &phi) const
        { return bi_(phi); }

        void init();

        const Vector2D &d() const
        { return d_; }

    protected:

        Vector2D d_;

        BilinearInterpolator bi_;
    };

    GhostCellImmersedBoundaryObject(const std::string &name,
                                    const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                    const std::shared_ptr<CellGroup> &solverCells);

    Type type() const
    { return GHOST_CELL; }

    void updateCells();

    Equation<Scalar> bcs(ScalarFiniteVolumeField &field) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField &field) const;

    Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    Equation<Scalar> pressureBcs(ScalarFiniteVolumeField &p) const;

    void computeForce(Scalar rho,
                      Scalar mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const Vector2D &g = Vector2D(0., 0.));

    void computeForce(const ScalarFiniteVolumeField &rho,
                      const ScalarFiniteVolumeField &mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const Vector2D &g = Vector2D(0., 0.));


    const std::vector<FixedStencil> &fixedStencils() const
    { return fixedStencils_; }

    const std::vector<ZeroGradientStencil> &zeroGradientStencils() const
    { return zeroGradientStencils_; }

protected:

    void constructStencils();

    std::vector<FixedStencil> fixedStencils_;

    std::vector<ZeroGradientStencil> zeroGradientStencils_;
};

#endif
