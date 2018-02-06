#include "QuadraticIbm.h"
#include "Matrix.h"
#include "QuadraticImmersedBoundaryObject.h"
#include "StaticMatrix.h"

Equation<Vector2D> qibm::div(const VectorFiniteVolumeField &phi,
                             VectorFiniteVolumeField &u,
                             const ImmersedBoundary &ib,
                             Scalar theta)
{
    Equation<Vector2D> eqn(u);
    std::vector<bool> isForcingCell(u.grid()->cells().size(), false);
    const VectorFiniteVolumeField &phi0 = phi.oldField(0);

    for (auto ibObj: ib)
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("qibm", "div", "immersed boundary objects must be of type \"quadratic\".");

        for (const Cell &cell: std::static_pointer_cast<QuadraticImmersedBoundaryObject>(ibObj)->forcingCells())
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar flux = dot(phi(nb.face()), nb.outwardNorm());
                Scalar flux0 = dot(phi0(nb.face()), nb.outwardNorm());

                eqn.add(cell, cell, theta * std::max(flux, 0.));
                eqn.addSource(cell, (1. - theta) * std::max(flux0, 0.) * u(cell));

                auto ibObj = ib.ibObj(nb.cell().centroid());

                if (ibObj)
                {
                    QuadraticIbmStencil st = QuadraticIbmStencil(nb, ib);
                    eqn.add(cell, st.cells(), std::valarray<Scalar>(std::min(flux, 0.) * st.coeffs()));
                    eqn.addSource(cell, std::min(flux, 0.) * st.src());
                }
                else
                {
                    eqn.add(cell, nb.cell(), theta * std::min(flux, 0.));
                    eqn.addSource(cell, (1. - theta) * std::min(flux0, 0.) * u(nb.cell()));
                }

                isForcingCell[cell.id()] = true;
            }
    }

    for (const Cell &cell: u.grid()->cellZone("fluid"))
    {
        if (isForcingCell[cell.id()])
            continue;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = theta * dot(phi(nb.face()), nb.outwardNorm());
            Scalar flux0 = (1. - theta) * dot(phi.oldField(0)(nb.face()), nb.outwardNorm());

            eqn.add(cell, cell, std::max(flux, 0.));
            eqn.add(cell, nb.cell(), std::min(flux, 0.));

            eqn.addSource(cell, std::max(flux0, 0.) * u(cell));
            eqn.addSource(cell, std::min(flux0, 0.) * u(nb.cell()));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = theta * dot(phi(bd.face()), bd.outwardNorm());
            Scalar flux0 = (1. - theta) * dot(phi(bd.face()), bd.outwardNorm());

            switch (u.boundaryType(bd.face()))
            {
                case VectorFiniteVolumeField::FIXED:
                    eqn.addSource(cell, flux * u(bd.face()));
                    eqn.addSource(cell, flux0 * u(bd.face()));
                    break;

                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                    eqn.add(cell, cell, flux);
                    eqn.addSource(cell, flux0 * u(cell));
                    break;

                case VectorFiniteVolumeField::SYMMETRY:
                    break;

                default:
                    throw Exception("qibm", "div", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}

Equation<Vector2D> qibm::laplacian(Scalar mu,
                                   VectorFiniteVolumeField &u,
                                   const ImmersedBoundary &ib,
                                   Scalar theta)
{
    Equation<Vector2D> eqn(u);
    std::vector<bool> isForcingCell(u.grid()->cells().size(), false);

    for (auto ibObj: ib)
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("qibm", "laplacian", "immersed boundary objects must be of type \"quadratic\".");

        for (const Cell &cell: std::static_pointer_cast<QuadraticImmersedBoundaryObject>(ibObj)->forcingCells())
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar flux = mu * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();

                eqn.add(cell, cell, -flux);

                auto ibObj = ib.ibObj(nb.cell().centroid());

                if (ibObj)
                {
                    QuadraticIbmStencil st = QuadraticIbmStencil(nb, ib);
                    eqn.add(cell, st.cells(), std::valarray<Scalar>(flux * st.coeffs()));
                    eqn.addSource(cell, flux * st.src());
                }
                else
                    eqn.add(cell, nb.cell(), flux);

                isForcingCell[cell.id()] = true;
            }
    }

    for (const Cell &cell: u.grid()->cellZone("fluid"))
    {
        if (isForcingCell[cell.id()])
            continue;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = mu * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
            eqn.add(cell, cell, -flux);
            eqn.add(cell, nb.cell(), flux);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = mu * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();

            switch (u.boundaryType(bd.face()))
            {
                case VectorFiniteVolumeField::FIXED:
                    eqn.add(cell, cell, -flux);
                    eqn.addSource(cell, flux * u(bd.face()));
                    break;

                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                case VectorFiniteVolumeField::SYMMETRY:
                    break;

                default:
                    throw Exception("qibm", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}

Equation<Vector2D> qibm::laplacian(const ScalarFiniteVolumeField &mu,
                                   VectorFiniteVolumeField &u,
                                   const ImmersedBoundary &ib,
                                   Scalar theta)
{
    Equation<Vector2D> eqn(u);
    std::vector<bool> isForcingCell(u.grid()->cells().size(), false);

    for (auto ibObj: ib)
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("qibm", "laplacian", "immersed boundary objects must be of type \"quadratic\".");

        for (const Cell &cell: std::static_pointer_cast<QuadraticImmersedBoundaryObject>(ibObj)->forcingCells())
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar flux = mu(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
                eqn.add(cell, cell, -flux);

                auto ibObj = ib.ibObj(nb.cell().centroid());

                if (ibObj)
                {
                    QuadraticIbmStencil st = QuadraticIbmStencil(nb, ib);
                    eqn.add(cell, st.cells(), std::valarray<Scalar>(flux * st.coeffs()));
                    eqn.addSource(cell, flux * st.src());
                }
                else
                    eqn.add(cell, nb.cell(), flux);

                isForcingCell[cell.id()] = true;
            }
    }

    for (const Cell &cell: u.grid()->cellZone("fluid"))
    {
        if (isForcingCell[cell.id()])
            continue;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = mu(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
            eqn.add(cell, cell, -flux);
            eqn.add(cell, nb.cell(), flux);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = mu(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();

            switch (u.boundaryType(bd.face()))
            {
                case VectorFiniteVolumeField::FIXED:
                    eqn.add(cell, cell, -flux);
                    eqn.addSource(cell, flux * u(bd.face()));
                    break;

                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                case VectorFiniteVolumeField::SYMMETRY:
                    break;

                default:
                    throw Exception("qibm", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return u;
}