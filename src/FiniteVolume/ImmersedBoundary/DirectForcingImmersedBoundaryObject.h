#ifndef DIRECT_FORCING_IMMERSED_BOUNDARY_OBJECT_H
#define DIRECT_FORCING_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"
#include "StaticMatrix.h"
#include "Matrix.h"

class DirectForcingImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    class Stencil
    {
    public:

        Stencil(const VectorFiniteVolumeField &u, const Cell &cell, const ImmersedBoundaryObject &ibObj);

        const Vector2D &uf() const
        { return uf_; }

        const std::vector<Scalar>& ipCoeffs() const
        { return ipCoeffs_; }

    protected:

        Stencil()
        {}

        Point2D bp_, ip_;

        Vector2D ub_, uip_, uf_;

        std::vector<Ref<const Cell>> ipCells_;

        std::vector<Scalar> ipCoeffs_;
    };

    class FieldExtensionStencil : public Stencil
    {
    public:
        FieldExtensionStencil(const VectorFiniteVolumeField &u, const Cell &cell, const ImmersedBoundaryObject &ibObj);
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
    {}

    virtual Equation<Vector2D> bcs(VectorFiniteVolumeField &field) const
    {}

    virtual Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const
    {}

    virtual void computeForce(Scalar rho,
                              Scalar mu,
                              const VectorFiniteVolumeField &u,
                              const ScalarFiniteVolumeField &p,
                              const Vector2D &g = Vector2D(0., 0.))
    {}

    void updateIbForce(const VectorFiniteVolumeField &u, Scalar timeStep, VectorFiniteVolumeField &fb);

private:

};


#endif
