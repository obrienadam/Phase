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

        Stencil(const Cell &cell,
                const ImmersedBoundaryObject &ibObj,
                const Vector2D& d);

        Stencil(const Cell& cell,
                const ImmersedBoundaryObject &ibObj);

        void initFixedBc();

        void initZeroGradientBc();

        const Cell &cell() const
        { return cell_; }

        const std::vector<Ref<const Cell>> &bpCells() const
        { return bpCells_; }

        const std::vector<Ref<const Cell>> &ipCells() const
        { return ipCells_; }

        const std::vector<Ref<const Cell>>& fixedBcCells() const
        { return fixedBcCells_; }

        const std::vector<Ref<const Cell>>& normalGradBcCells() const
        { return zeroGradBcCells_; }

        const std::vector<Scalar> &fixedBcCoeffs() const
        { return fixedBcCoeffs_; }

        const std::vector<Scalar> &normalGradBcCoeffs() const
        { return zeroGradBcCoeffs_; }

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

        StaticMatrix<4, 4> Aip_, Abp_;

        Ref<const Cell> cell_;

        std::vector<Ref<const Cell>> ipCells_, bpCells_, fixedBcCells_, zeroGradBcCells_;

        std::vector<Scalar> fixedBcCoeffs_, zeroGradBcCoeffs_;

        Point2D ip_, bp_, nw_, d_;
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


    const std::vector<Stencil> &stencils() const
    { return stencils_; }

protected:

    void constructStencils();

    std::vector<Stencil> stencils_;
};

#endif
