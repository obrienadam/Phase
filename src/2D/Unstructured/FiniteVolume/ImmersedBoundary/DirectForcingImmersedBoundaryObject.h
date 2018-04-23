#ifndef PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_OBJECT_H
#define PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_OBJECT_H

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

    class QuadraticLsStencil
    {
    public:

        QuadraticLsStencil(const VectorFiniteVolumeField &u,
                           const Cell &cell,
                           const ImmersedBoundaryObject &ibObj);

        const Point2D &bp() const
        {return bp_; }

        const Vector2D &uf() const
        { return uf_; }

    protected:

        Point2D bp_;

        Vector2D ub_, uf_;

    };

    //- Constructors, one for circles, another for polygons
    DirectForcingImmersedBoundaryObject(const std::string &name,
                                        const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                        const std::shared_ptr<CellGroup> &solverCells);

    Type type() const
    { return DIRECT_FORCING; }

    void updateCells();

    virtual FiniteVolumeEquation<Scalar> bcs(ScalarFiniteVolumeField &field) const
    { return FiniteVolumeEquation<Scalar>(field); }

    virtual FiniteVolumeEquation<Vector2D> bcs(VectorFiniteVolumeField &field) const
    { return FiniteVolumeEquation<Vector2D>(field); }

    virtual FiniteVolumeEquation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

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
