#include "AxisymmetricLaplacian.h"

Equation<Scalar> axi::laplacian(Scalar gamma, ScalarFiniteVolumeField &phi, Scalar theta)
{
    Equation<Scalar> eqn(phi);

    for (const Cell &cell: phi.cells())
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.face().polarOutwardNorm(cell.centroid());
            Scalar flux = gamma * dot(nb.rCellVec(), sf) / nb.rCellVec().magSqr();
            eqn.add(cell, cell, -flux);
            eqn.add(cell, nb.cell(), flux);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.face().polarOutwardNorm(cell.centroid());
            Scalar flux = gamma * dot(bd.rFaceVec(), sf) / bd.rFaceVec().magSqr();

            switch (phi.boundaryType(bd.face()))
            {
                case ScalarFiniteVolumeField::FIXED:
                    eqn.add(cell, cell, -flux);
                    eqn.addSource(cell, flux * phi(bd.face()));
                    break;

                case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                case ScalarFiniteVolumeField::SYMMETRY:
                    break;

                default:
                    throw Exception("axi", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}

Equation<Vector2D> axi::vectorLaplacian(Scalar gamma, VectorFiniteVolumeField &phi, Scalar theta)
{
    Equation<Vector2D> eqn(phi);

    for (const Cell &cell: phi.cells())
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.face().polarOutwardNorm(cell.centroid());
            Scalar flux = gamma * dot(nb.rCellVec(), sf) / nb.rCellVec().magSqr();

            eqn.add(cell, cell, -theta * flux);
            eqn.add(cell, nb.cell(), theta * flux);
            eqn.addSource(cell, (1. - theta) * flux * (phi(nb.cell()) - phi(cell)));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.face().polarOutwardNorm(cell.centroid());
            Scalar flux = gamma * dot(bd.rFaceVec(), sf) / bd.rFaceVec().magSqr();

            switch (phi.boundaryType(bd.face()))
            {
                case VectorFiniteVolumeField::FIXED:
                    eqn.add(cell, cell, -theta * flux);
                    eqn.addSource(cell, theta * flux * phi(bd.face()));
                    eqn.addSource(cell, (1. - theta) * flux * (phi(bd.face()) - phi(cell)));
                    break;

                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                case VectorFiniteVolumeField::SYMMETRY:
                    break;

                default:
                    throw Exception("axi", "vectorLaplacian", "unrecognized or unspecified boundary type.");
            }
        }

        Scalar vol = cell.polarVolume();

        eqn.add(cell, cell, Vector2D(-theta * gamma * vol / std::pow(cell.centroid().x, 2), 0.));
        eqn.addSource(cell, Vector2D(
                -(1. - theta) * gamma * vol * phi(cell).x / std::pow(cell.centroid().x, 2), 0.));
    }

    return eqn;
}