#include "Slic.h"

Equation<Scalar> slic::div(const VectorFiniteVolumeField &u,
                           const VectorFiniteVolumeField &gradGamma,
                           ScalarFiniteVolumeField &gamma,
                           Scalar timeStep)
{
    Equation<Scalar> eqn(gamma);

    for (const Cell &cell: gamma.grid().cellZone("fluid"))
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D m = -gradGamma(cell).unitVec();


        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = dot(u(bd.face()), bd.outwardNorm());
            switch (gamma.boundaryType(bd.face()))
            {
                case ScalarFiniteVolumeField::FIXED:
                    eqn.addSource(cell, flux * gamma(bd.face()));
                    break;

                case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                    eqn.add(cell, cell, flux);
                    eqn.addSource(cell, flux*gamma(bd.face()));
                    break;

                case ScalarFiniteVolumeField::SYMMETRY:
                    break;

                default:
                    throw Exception("slic", "div", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}