#ifndef GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H
#define GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"
#include "StaticMatrix.h"

class GhostCellImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    class Stencil
    {
    public:

        Stencil(const Cell &cell) : cell_(std::cref(cell))
        {}

        Stencil(const Cell &cell,
                const ImmersedBoundaryObject &ibObj);

        Stencil(const Cell &cell,
                const ImmersedBoundaryObject &ibObj,
                const Vector2D &r);

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

        Scalar bpValue(const ScalarFiniteVolumeField &phi) const;

        Vector2D bpValue(const VectorFiniteVolumeField &u) const;

        Vector2D bpGrad(const ScalarFiniteVolumeField &phi) const;

        Tensor2D bpGrad(const VectorFiniteVolumeField &u) const;

    protected:

        void initMatrices();

        StaticMatrix<4, 4> Aip_, Abp_;

        Ref<const Cell> cell_;

        std::vector<Ref<const Cell>> ipCells_, bpCells_, cells_;

        std::vector<Scalar> coeffs_;

        Point2D ip_, bp_, nw_;
    };

    class FixedStencil : public Stencil
    {
    public:

        FixedStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj)
                :
                Stencil(cell, ibObj)
        {
            init();
        }

        void init();
    };

    class ZeroGradientStencil : public Stencil
    {
    public:

        ZeroGradientStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj)
                :
                Stencil(cell, ibObj)
        {
            init();
        }

        ZeroGradientStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj, const Vector2D &r)
                :
                Stencil(cell, ibObj, r)
        {
            init();
        }

        void init();
    };

    class ForcingStencil : public Stencil
    {
    public:

        ForcingStencil(const Cell &cell,
                       const ImmersedBoundaryObject &ibObj);

        void init();

    private:
    };

    GhostCellImmersedBoundaryObject(const std::string &name,
                                    Label id,
                                    const ImmersedBoundary &ib,
                                    const std::shared_ptr<FiniteVolumeGrid2D> &grid);

    Type type() const
    { return GHOST_CELL; }

    void updateCells();

    Equation<Scalar> bcs(ScalarFiniteVolumeField &field) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField &field) const;

    Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

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


    const std::vector<std::shared_ptr<Stencil>> &fixedStencils() const
    { return fixedStencils_; }

    const std::vector<ZeroGradientStencil> &zeroGradientStencils() const
    { return zeroGradientStencils_; }

protected:

    void constructStencils();

    std::vector<std::shared_ptr<Stencil>> fixedStencils_;

    std::vector<ZeroGradientStencil> zeroGradientStencils_;
};

#endif
