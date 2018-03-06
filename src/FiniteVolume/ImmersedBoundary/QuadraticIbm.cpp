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

    for (const Cell &cell: u.grid()->cellZone("fluid"))
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = dot(phi(nb.face()), nb.outwardNorm());

            eqn.add(cell, cell, std::max(flux, 0.));

            auto lIbObj = ib.ibObj(cell.centroid());
            auto rIbObj = ib.ibObj(nb.cell().centroid());
            auto ibObj = lIbObj ? lIbObj : rIbObj;

            if (ibObj)
            {
                QuadraticIbmStencil st = QuadraticIbmStencil(nb, ib);
                eqn.add(cell, st.cells(), std::valarray<Scalar>(std::min(flux, 0.) * st.coeffs()));
                eqn.addSource(cell, std::min(flux, 0.) * st.src());
            }
            else
            {
                eqn.add(cell, nb.cell(), theta * std::min(flux, 0.));
            }
        }
        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = dot(phi(bd.face()), bd.outwardNorm());

            switch (u.boundaryType(bd.face()))
            {
                case VectorFiniteVolumeField::FIXED:eqn.addSource(cell, flux * u(bd.face()));
                    break;

                case VectorFiniteVolumeField::NORMAL_GRADIENT:eqn.add(cell, cell, flux);
                    break;

                case VectorFiniteVolumeField::SYMMETRY:break;

                default:throw Exception("qibm", "div", "unrecognized or unspecified boundary type.");
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

    for (const Cell &cell: u.grid()->cellZone("fluid"))
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = mu * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();

            eqn.add(cell, cell, -flux);

            auto lIbObj = ib.ibObj(cell.centroid());
            auto rIbObj = ib.ibObj(nb.cell().centroid());
            auto ibObj = lIbObj ? lIbObj : rIbObj;

            if (ibObj)
            {
                QuadraticIbmStencil st = QuadraticIbmStencil(nb, ib);
                eqn.add(cell, st.cells(), std::valarray<Scalar>(flux * st.coeffs()));
                eqn.addSource(cell, flux * st.src());
            }
            else
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
                case VectorFiniteVolumeField::SYMMETRY:break;

                default:throw Exception("qibm", "laplacian", "unrecognized or unspecified boundary type.");
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

    for (const Cell &cell: u.grid()->cellZone("fluid"))
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = mu(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();

            eqn.add(cell, cell, -flux);

            auto lIbObj = ib.ibObj(cell.centroid());
            auto rIbObj = ib.ibObj(nb.cell().centroid());
            auto ibObj = lIbObj ? lIbObj : rIbObj;

            if (ibObj)
            {
                QuadraticIbmStencil st = QuadraticIbmStencil(nb, ib);
                eqn.add(cell, st.cells(), std::valarray<Scalar>(flux * st.coeffs()));
                eqn.addSource(cell, flux * st.src());
            }
            else
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
                case VectorFiniteVolumeField::SYMMETRY:break;

                default:throw Exception("qibm", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}

void qibm::interpolateFaces(VectorFiniteVolumeField &u, const ImmersedBoundary &ib)
{
    for (const Face &face: u.grid()->interiorFaces())
    {
        auto lIbObj = ib.ibObj(face.lCell().centroid());
        auto rIbObj = ib.ibObj(face.rCell().centroid());

        if (lIbObj == rIbObj)
        {
            Scalar g = face.volumeWeight();
            u(face) = g * u(face.lCell()) + (1. - g) * u(face.rCell());
            continue;
        }
        else
        {
            u(face) = Vector2D(0., 0.);
            continue;
        }

        auto ibObj = lIbObj ? lIbObj : rIbObj;

        Vector2D eta = (face.rCell().centroid() - face.lCell().centroid()).unitVec();
        Point2D xb = ibObj->intersectionLine(face.lCell().centroid(), face.rCell().centroid()).ptB();

        Scalar e[] = {
                dot(face.lCell().centroid(), eta),
                dot(xb, eta),
                dot(face.rCell().centroid(), eta)
        };

        Vector2D v[] = {
                u(face.lCell()),
                ibObj->velocity(xb),
                u(face.rCell())
        };

        auto A = StaticMatrix<3, 3>(
                {
                        e[0] * e[0], e[0], 1.,
                        e[1] * e[1], e[1], 1.,
                        e[2] * e[2], e[2], 1.
                });

        auto b = StaticMatrix<3, 2>(
                {
                        v[0].x, v[0].y,
                        v[1].x, v[1].y,
                        v[2].x, v[2].y,
                });

        auto c = solve(A, b);

        Scalar etaf = dot(face.centroid(), eta);

        Vector2D uf = Vector2D(
                c(0, 0) * etaf * etaf + c(1, 0) * etaf + c(2, 0),
                c(0, 1) * etaf * etaf + c(1, 1) * etaf + c(2, 1)
        );

        std::cout << uf << std::endl;
        u(face) = uf;
    }

    u.interpolateFaces();
}