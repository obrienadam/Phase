#include "AxisymmetricLaplacian.h"

FiniteVolumeEquation<Scalar> axi::laplacian(Scalar gamma,
                                            ScalarFiniteVolumeField &phi)
{
    FiniteVolumeEquation<Scalar> eqn(phi);

    for (const Cell &cell: phi.cells())
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.polarOutwardNorm();
            Scalar flux = gamma * dot(nb.rc(), sf) / nb.rc().magSqr();
            eqn.add(cell, cell, -flux);
            eqn.add(cell, nb.cell(), flux);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.polarOutwardNorm();
            Scalar flux = gamma * dot(bd.rf(), sf) / bd.rf().magSqr();

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

FiniteVolumeEquation<Scalar> axi::laplacian(const ScalarFiniteVolumeField &gamma, ScalarFiniteVolumeField &phi)
{
    FiniteVolumeEquation<Scalar> eqn(phi, 5);

    for (const Cell &cell: phi.cells())
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.polarOutwardNorm();
            Scalar flux = gamma(nb.face()) * dot(nb.rc(), sf) / nb.rc().magSqr();
            eqn.add(cell, cell, -flux);
            eqn.add(cell, nb.cell(), flux);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.polarOutwardNorm();
            Scalar flux = gamma(bd.face()) * dot(bd.rf(), sf) / bd.rf().magSqr();

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

FiniteVolumeEquation<Vector2D> axi::laplacian(Scalar gamma,
                                              VectorFiniteVolumeField &u,
                                              Scalar theta)
{
    FiniteVolumeEquation<Vector2D> eqn(u);
    const VectorFiniteVolumeField &u0 = u.oldField(0);

    for (const Cell &cell: u.cells())
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.polarOutwardNorm();
            Scalar flux = gamma * dot(nb.rc(), sf) / nb.rc().magSqr();

            eqn.add(cell, cell, -theta * flux);
            eqn.add(cell, nb.cell(), theta * flux);
            eqn.addSource(cell, (1. - theta) * flux * (u0(nb.cell()) - u0(cell)));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.polarOutwardNorm();

            Scalar flux = gamma * dot(bd.rf(), sf) / bd.rf().magSqr();

            switch (u.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.add(cell, cell, -theta * flux);
                eqn.addSource(cell, theta * flux * u0(bd.face()));
                eqn.addSource(cell, (1. - theta) * flux * (u0(bd.face()) - u0(cell)));
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                break;
            case VectorFiniteVolumeField::SYMMETRY:
            {
                Vector2D tw = bd.outwardNorm().tangentVec().unitVec();

                eqn.add(cell, cell, -theta * flux);
                eqn.add(cell, cell, theta * flux * outer(tw, tw));
                eqn.addSource(cell, (1. - theta) * flux * (dot(u0(cell), tw) * tw - u0(cell)));
            }
                break;

            default:
                throw Exception("axi", "vectorLaplacian", "unrecognized or unspecified boundary type.");
            }
        }

        //- Non-polar volume since we just want the rz face
        //Scalar r = cell.centroid().x;
        //eqn.add(cell, cell, -gamma * cell.volume() / r * Vector2D(theta, 0.));
        //eqn.addSource(cell, -gamma * cell.volume() * u0(cell).x / r * Vector2D(1. - theta, 0.));
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> axi::laplacian(const ScalarFiniteVolumeField &gamma,
                                              VectorFiniteVolumeField &u,
                                              Scalar theta)
{
    FiniteVolumeEquation<Vector2D> eqn(u);
    const ScalarFiniteVolumeField &gamma0 = gamma.oldField(0);
    const VectorFiniteVolumeField &u0 = u.oldField(0);

    for (const Cell &cell: u.cells())
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.polarOutwardNorm();
            Scalar flux = gamma(nb.face()) * dot(nb.rc(), sf) / nb.rc().magSqr();
            Scalar flux0 = gamma0(nb.face()) * dot(nb.rc(), sf) / nb.rc().magSqr();

            eqn.add(cell, cell, -theta * flux);
            eqn.add(cell, nb.cell(), theta * flux);
            eqn.addSource(cell, (1. - theta) * flux0 * (u0(nb.cell()) - u0(cell)));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.polarOutwardNorm();

            Scalar flux = gamma(bd.face()) * dot(bd.rf(), sf) / bd.rf().magSqr();
            Scalar flux0 = gamma0(bd.face()) * dot(bd.rf(), sf) / bd.rf().magSqr();

            switch (u.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.add(cell, cell, -theta * flux);
                eqn.addSource(cell, theta * flux * u0(bd.face()));
                eqn.addSource(cell, (1. - theta) * flux0 * (u0(bd.face()) - u0(cell)));
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                break;
            case VectorFiniteVolumeField::SYMMETRY:
            {
                Vector2D tw = bd.outwardNorm().tangentVec().unitVec();

                eqn.add(cell, cell, -theta * flux);
                eqn.add(cell, cell, theta * flux * outer(tw, tw));
                eqn.addSource(cell, (1. - theta) * flux0 * (dot(u0(cell), tw) * tw - u0(cell)));
            }
                break;

            default:
                throw Exception("axi", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }

        //- Non-polar volume since we just want the rz face
//        Scalar r = cell.centroid().x;
//        eqn.add(cell, cell, -gamma(cell) * cell.volume() / r * Vector2D(theta, 0.));
//        eqn.addSource(cell, -gamma0(cell) * cell.volume() * u0(cell).x / r * Vector2D(1. - theta, 0.));
    }

    return eqn;
}
