#include "QuadraticIbm.h"
#include "Matrix.h"
#include "QuadraticImmersedBoundaryObject.h"

Equation<Vector2D> qibm::div(const VectorFiniteVolumeField &phi,
                             VectorFiniteVolumeField &u,
                             const ImmersedBoundary &ib,
                             Scalar theta)
{
    Equation<Vector2D> eqn(u);
    std::vector<bool> isForcingCell(u.grid().cells().size(), false);

    for (auto ibObj: ib)
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("qibm", "div", "immersed boundary objects must be of type \"quadratic\".");

        for (const Cell &cell: std::static_pointer_cast<QuadraticImmersedBoundaryObject>(ibObj)->forcingCells())
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar flux = dot(phi(nb.face()), nb.outwardNorm());
                eqn.add(cell, cell, std::max(flux, 0.));

                auto ibObj = ib.ibObj(nb.cell().centroid());

                if (ibObj)
                {
                    QuadraticIbmStencil st = QuadraticIbmStencil(cell, nb.cell(), ib, std::min(flux, 0.));
                    eqn.add(cell, st.cells(), st.coeffs());
                    eqn.addSource(cell, st.src());
                }
                else
                    eqn.add(cell, nb.cell(), std::min(flux, 0.));

                isForcingCell[cell.id()] = true;
            }
    }

    for (const Cell &cell: u.grid().cellZone("fluid"))
    {
        if(isForcingCell[cell.id()])
            continue;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = dot(phi(nb.face()), nb.outwardNorm());
            eqn.add(cell, cell, std::max(flux, 0.));
            eqn.add(cell, nb.cell(), std::min(flux, 0.));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = dot(phi(bd.face()), bd.outwardNorm());

            switch (u.boundaryType(bd.face()))
            {
                case VectorFiniteVolumeField::FIXED:
                    eqn.addSource(cell, flux * u(bd.face()));
                    break;

                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                    eqn.add(cell, cell, flux);
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
                                   const ImmersedBoundary &ib)
{
    Equation<Vector2D> eqn(u);
    std::vector<bool> isForcingCell(u.grid().cells().size(), false);

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
                    QuadraticIbmStencil st = QuadraticIbmStencil(cell, nb.cell(), ib, flux);
                    eqn.add(cell, st.cells(), st.coeffs());
                    eqn.addSource(cell, st.src());
                }
                else
                    eqn.add(cell, nb.cell(), flux);

                isForcingCell[cell.id()] = true;
            }
    }

    for (const Cell &cell: u.grid().cellZone("fluid"))
    {
        if(isForcingCell[cell.id()])
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
                                   const ImmersedBoundary &ib)
{
    Equation<Vector2D> eqn(u);
    std::vector<bool> isForcingCell(u.grid().cells().size(), false);

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
                    QuadraticIbmStencil st = QuadraticIbmStencil(cell, nb.cell(), ib, flux);
                    eqn.add(cell, st.cells(), st.coeffs());
                    eqn.addSource(cell, st.src());
                }
                else
                    eqn.add(cell, nb.cell(), flux);

                isForcingCell[cell.id()] = true;
            }
    }

    for (const Cell &cell: u.grid().cellZone("fluid"))
    {
        if(isForcingCell[cell.id()])
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