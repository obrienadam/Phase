#ifndef DIRECT_FORCING_IMMERSED_BOUNDARY_OBJECT_H
#define DIRECT_FORCING_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"

class DirectForcingImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    class Stencil
    {
    public:

        Stencil(const VectorFiniteVolumeField &u,
                const Cell &cell,
                const ImmersedBoundaryObject &ibObj);

        const Point2D &bp() const
        { return bp_; }

        const Point2D &ip() const
        { return ip_; }

        const Vector2D &uf() const
        { return uf_; }

    protected:

        Stencil()
        {}

        Point2D bp_, ip_;

        Vector2D ub_, uip_, uf_;
    };

    class FieldExtensionStencil : public Stencil
    {
    public:
        FieldExtensionStencil(const VectorFiniteVolumeField &u,
                              const Cell &cell,
                              const ImmersedBoundaryObject &ibObj);

        const std::vector<Ref<const Cell>> &cells() const
        { return cells_; }

        const std::vector<Scalar> &coeffs() const
        { return coeffs_; }

    private:

        std::vector<Ref<const Cell>> cells_;

        std::vector<Scalar> coeffs_;
    };

    class PressureFieldExtensionStencil
    {
    public:
        PressureFieldExtensionStencil(const ScalarFiniteVolumeField &p,
                                      const Cell &cell,
                                      const ImmersedBoundaryObject &ibObj);

        const Point2D &bp() const
        { return bp_; }

        const Vector2D &n() const
        { return n_; }

        const std::vector<Ref<const Cell>> &cells() const
        { return cells_; }

        const std::vector<Scalar> &coeffs() const
        { return coeffs_; }

    protected:

        Point2D bp_, ip_;

        Vector2D n_;

        std::vector<Ref<const Cell>> cells_;

        std::vector<Scalar> coeffs_;
    };

    class QuadraticStencil
    {
    public:

        QuadraticStencil(const VectorFiniteVolumeField &u, const Cell &cell, const ImmersedBoundaryObject &ibObj);

        const Vector2D &uf() const
        { return uf_; }

    private:

        Vector2D uf_;
    };

    //- Constructors, one for circles, another for polygons
    DirectForcingImmersedBoundaryObject(const std::string &name,
                                        Label id,
                                        const ImmersedBoundary &ib,
                                        const std::shared_ptr<FiniteVolumeGrid2D> &grid);

    Type type() const
    { return DIRECT_FORCING; }

    void updateCells();

    virtual Equation<Scalar> bcs(ScalarFiniteVolumeField &field) const
    { return Equation<Scalar>(field); }

    virtual Equation<Vector2D> bcs(VectorFiniteVolumeField &field) const
    { return Equation<Vector2D>(field); }

    virtual Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const
    { return Equation<Vector2D>(u); }

    virtual void computeForce(Scalar rho,
                              Scalar mu,
                              const VectorFiniteVolumeField &u,
                              const ScalarFiniteVolumeField &p,
                              const Vector2D &g = Vector2D(0., 0.))
    {}

    void computeBoundaryForcing(const VectorFiniteVolumeField &u,
                                Scalar timeStep,
                                VectorFiniteVolumeField &fb) const;

    const CellGroup &forcingCells() const
    { return forcingCells_; }

    const CellGroup &pseudoForcingCells() const
    { return pseudoFluidCells_; }

private:

    CellGroup forcingCells_, pseudoFluidCells_;

};


#endif
