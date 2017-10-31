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
    std::vector<bool> isForcingCell(u.grid().cells().size(), false);
    const VectorFiniteVolumeField &phi0 = phi.oldField(0);

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
        if (isForcingCell[cell.id()])
            continue;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = theta * dot(phi(nb.face()), nb.outwardNorm());
            Scalar flux0 = (1. - theta) * dot(phi.oldField(0)(nb.face()), nb.outwardNorm());

            eqn.add(cell, cell, std::max(flux, 0.));
            eqn.add(cell, nb.cell(), std::min(flux, 0.));

            eqn.addSource(cell, std::max(flux0, 0.) * u.oldField(0)(cell));
            eqn.addSource(cell, std::min(flux0, 0.) * u.oldField(0)(nb.cell()));
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

Equation<Scalar>
qibm::pressureEquation(Scalar rho, Scalar timeStep, ScalarFiniteVolumeField &p, const ImmersedBoundary &ib)
{
    Equation<Scalar> eqn(p);
    std::vector<bool> isForcingCell(p.grid().cells().size(), false);

    for (auto ibObj: ib)
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("qibm", "pressureEquation", "immersed boundary objects must be of type \"quadratic\".");

        for (const Cell &cell: std::static_pointer_cast<QuadraticImmersedBoundaryObject>(ibObj)->forcingCells())
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar flux = timeStep / rho * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
                eqn.add(cell, cell, -flux);

                auto ibObj = ib.ibObj(nb.cell().centroid());

                if (ibObj)
                {
                    const Cell &cell0 = nb.cell();
                    const Cell &cell1 = cell;
                    const Cell &cell2 = p.grid().globalActiveCells().nearestItem(
                            2 * cell1.centroid() - cell0.centroid());

                    LineSegment2D ln = ibObj->intersectionLine(cell1.centroid(), cell0.centroid());
                    Vector2D eta = (cell0.centroid() - cell2.centroid());

                    Scalar eta0 = dot(cell0.centroid(), eta);
                    Scalar eta1 = dot(cell1.centroid(), eta);
                    Scalar eta2 = dot(cell2.centroid(), eta);
                    Scalar etaB = dot(ln.ptB(), eta);

                    auto A = inverse(StaticMatrix<3, 3>(
                            {
                                    eta1 * eta1, eta1, 1.,
                                    eta2 * eta2, eta2, 1.,
                                    2. * etaB, 1., 0.
                            }));

                    auto c = StaticMatrix<1, 3>({eta0 * eta0, eta0, 1.}) * A;

                    eqn.add(cell, cell1, c(0, 0) * flux);
                    eqn.add(cell, cell2, c(0, 1) * flux);

                    Scalar DuDt = -1. / rho * dot(ibObj->acceleration(ln.ptB()), eta);
                    //eqn.addSource(cell, c(0, 2) * flux * DuDt);
                }
                else
                    eqn.add(cell, nb.cell(), flux);

                isForcingCell[cell.id()] = true;
            }
    }

    for (const Cell &cell: p.grid().localActiveCells())
    {
        if (isForcingCell[cell.id()])
            continue;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = timeStep / rho * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
            eqn.add(cell, cell, -flux);
            eqn.add(cell, nb.cell(), flux);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = timeStep / rho * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();

            switch (p.boundaryType(bd.face()))
            {
                case VectorFiniteVolumeField::FIXED:
                    eqn.add(cell, cell, -flux);
                    eqn.addSource(cell, flux * p(bd.face()));
                    break;

                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                case VectorFiniteVolumeField::SYMMETRY:
                    break;

                default:
                    throw Exception("qibm", "pressureEquation", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}

void qibm::computeFaceVelocities(VectorFiniteVolumeField &u, const ImmersedBoundary &ib)
{
    for (const Face &face: u.grid().interiorFaces())
    {
        auto fIbObj = ib.ibObj(face.centroid());

        if (fIbObj)
        {
            u(face) = fIbObj->velocity(face.centroid());
            continue;
        }

        auto lIbObj = ib.ibObj(face.lCell().centroid());
        auto rIbObj = ib.ibObj(face.rCell().centroid());
        Scalar g = face.distanceWeight();
        Vector2D eta = (face.rCell().centroid() - face.lCell().centroid()).unitVec();

        if ((bool) lIbObj == (bool) rIbObj) //- Either fluid cell or completely in solid
        {
            u(face) = g * u(face.lCell()) + (1. - g) * u(face.rCell());
            continue;
        }

        LineSegment2D ln = lIbObj ? lIbObj->intersectionLine(face.rCell().centroid(), face.lCell().centroid())
                                  : rIbObj->intersectionLine(face.lCell().centroid(), face.rCell().centroid());

        const Cell &cell0 = lIbObj ? face.lCell() : face.rCell();
        const Cell &cell1 = lIbObj ? face.rCell() : face.lCell();
        const Cell &cell2 = u.grid().globalActiveCells().nearestItem(2 * cell1.centroid() - cell0.centroid());

        Scalar eta1 = dot(cell1.centroid(), eta);
        Scalar eta2 = dot(cell2.centroid(), eta);
        Scalar etaB = dot(ln.ptB(), eta);
        Vector2D vb = lIbObj ? lIbObj->velocity(ln.ptB()) : rIbObj->velocity(ln.ptB());

        auto c = solve(StaticMatrix<3, 3>({
                                                  eta1 * eta1, eta1, 1.,
                                                  eta2 * eta2, eta2, 1.,
                                                  etaB * etaB, etaB, 1.
                                          }),
                       StaticMatrix<3, 2>({
                                                  u(cell1).x, u(cell1).y,
                                                  u(cell2).x, u(cell2).y,
                                                  vb.x, vb.y
                                          }));

        Scalar etaF = dot(face.centroid(), eta);

        u(face) = Vector2D(
                c(0, 0) * etaF * etaF + c(1, 0) * etaF + c(2, 0),
                c(0, 1) * etaF * etaF + c(1, 1) * etaF + c(2, 1)
        );
    }

    u.setBoundaryFaces();
}