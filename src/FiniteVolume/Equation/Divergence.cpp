#include "Divergence.h"

Equation<Vector2D> fv::div(const VectorFiniteVolumeField &phiU,
                           const JacobianField &gradU,
                           VectorFiniteVolumeField &u)
{
    Equation<Vector2D> eqn(u);

    for (const Cell &cell: u.grid()->cellZone("fluid"))
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = dot(phiU(nb.face()), nb.outwardNorm());

            if (flux > 0.)
            {
                eqn.add(cell, cell, flux);
                eqn.addSource(cell, flux * dot(gradU(cell), nb.rFaceVec()));
            }
            else
            {
                eqn.add(cell, nb.cell(), flux);
                eqn.addSource(cell, flux * dot(gradU(nb.cell()), nb.face().centroid() - nb.cell().centroid()));
            }
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = dot(phiU(bd.face()), bd.outwardNorm());

            switch (u.boundaryType(bd.face()))
            {
                case VectorFiniteVolumeField::FIXED:
                    eqn.addSource(cell, flux * u(bd.face()));
                    break;

                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                    eqn.add(cell, cell, flux);
                    eqn.addSource(cell, flux * dot(gradU(cell), bd.rFaceVec()));
                    break;

                case VectorFiniteVolumeField::SYMMETRY:
                    break;

                default:
                    throw Exception("fv", "div", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}