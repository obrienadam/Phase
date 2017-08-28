#include "QuadraticIbm.h"
#include "Matrix.h"
#include "Tensor2D.h"

Equation<Vector2D> qibm::div(const VectorFiniteVolumeField &phi,
                             VectorFiniteVolumeField &u,
                             const ImmersedBoundary& ib)
{
    Equation<Vector2D> eqn(u);
    const CellZone& fluid = u.grid().cellZone("fluid");

    for(const Cell& cell: fluid)
    {
        for(const InteriorLink& nb: cell.neighbours())
        {
            Scalar flux = dot(phi(nb.face()), nb.outwardNorm());
            eqn.add(cell, cell, std::max(flux, 0.));

            auto ibObj = ib.ibObj(nb.cell().centroid());

            if(ibObj)
            {
                const Cell &stCell = fluid.nearestItems(2*cell.centroid() - nb.cell().centroid(), 1)[0];
                LineSegment2D ln = ibObj->intersectionLine(LineSegment2D(cell.centroid(), nb.cell().centroid()));

                Vector2D eta = nb.rCellVec().unitVec();
                Scalar eta1 = dot(stCell.centroid(), eta);
                Scalar eta2 = dot(cell.centroid(), eta);
                Scalar eta3 = dot(ln.ptB(), eta);
                Scalar eta4 = dot(nb.cell().centroid(), eta);

                Matrix c = Matrix(1, 3, {eta4*eta4, eta4, 1.})*inverse(Matrix(3, 3, {
                        eta1*eta1, eta1, 1.,
                        eta2*eta2, eta2, 1.,
                        eta3*eta3, eta3, 1.
                }));

                eqn.add(cell, stCell, c(0, 0)*std::min(flux, 0.));
                eqn.add(cell, cell, c(0, 1)*std::min(flux, 0.));
                eqn.addSource(cell, c(0, 2)*std::min(flux, 0.)*ibObj->velocity(ln.ptB()));
            }
            else
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
                                   const ImmersedBoundary& ib)
{
    Equation<Vector2D> eqn(u);
    const CellZone& fluid = u.grid().cellZone("fluid");

    for(const Cell& cell: fluid)
    {
        for(const InteriorLink& nb: cell.neighbours())
        {
            Scalar flux = mu * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
            eqn.add(cell, cell, -flux);

            auto ibObj = ib.ibObj(nb.cell().centroid());

            if(ibObj)
            {
                const Cell &stCell = fluid.nearestItems(2*cell.centroid() - nb.cell().centroid(), 1)[0];
                LineSegment2D ln = ibObj->intersectionLine(LineSegment2D(cell.centroid(), nb.cell().centroid()));

                Vector2D eta = nb.rCellVec().unitVec();
                Scalar eta1 = dot(stCell.centroid(), eta);
                Scalar eta2 = dot(cell.centroid(), eta);
                Scalar eta3 = dot(ln.ptB(), eta);
                Scalar eta4 = dot(nb.cell().centroid(), eta);

                Matrix c = Matrix(1, 3, {eta4*eta4, eta4, 1.})*inverse(Matrix(3, 3, {
                        eta1*eta1, eta1, 1.,
                        eta2*eta2, eta2, 1.,
                        eta3*eta3, eta3, 1.
                }));

                eqn.add(cell, stCell, c(0, 0)*flux);
                eqn.add(cell, cell, c(0, 1)*flux);
                eqn.addSource(cell, c(0, 2)*flux*ibObj->velocity(ln.ptB()));
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
                case VectorFiniteVolumeField::SYMMETRY:
                    break;

                default:
                    throw Exception("qibm", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}

void qibm::computeFaceVelocities(VectorFiniteVolumeField& u, const ImmersedBoundary &ib)
{
    const CellZone& fluid = u.grid().cellZone("fluid");

    auto ibInterpolate = [&u, &fluid](const Cell& cell2, const Cell& cell3, const ImmersedBoundaryObject& ibObj) {
        const Cell& cell1 = fluid.nearestItems(2*cell2.centroid() - cell3.centroid(), 1)[0];
        LineSegment2D ln = ibObj.intersectionLine(LineSegment2D(cell2.centroid(), cell3.centroid()));

        Vector2D eta = (cell3.centroid() - cell2.centroid()).unitVec();
        Scalar eta1 = dot(cell1.centroid(), eta);
        Scalar eta2 = dot(cell2.centroid(), eta);
        Scalar eta3 = dot(ln.ptB(), eta);

        Matrix A(3, 3, {
                eta1*eta1, eta1, 1.,
                eta2*eta2, eta2, 1.,
                eta3*eta3, eta3, 1.
        });

        Vector2D u1 = u(cell1);
        Vector2D u2 = u(cell2);
        Vector2D u3 = ibObj.velocity(ln.ptB());

        Matrix bu = Matrix(3, 1, {u1.x, u2.x, u3.x});
        Matrix bv = Matrix(3, 1, {u1.y, u2.y, u3.y});

        Matrix cu = solve(A, bu);
        Matrix cv = solve(A, bv);

        Scalar eta4 = dot(cell3.centroid(), eta);
        Matrix x = Matrix(1, 3, {eta4*eta4, eta4, 1.});

        Vector2D u4 = Vector2D(
                (x*cu)(0, 0),
                (x*cv)(0, 0)
        );

        Scalar g = cell3.volume()/(cell2.volume() + cell3.volume());

        return g*u(cell2) + (1. - g)*u4;
    };

    for(const Face& f: u.grid().interiorFaces())
    {
        auto lCellIbObj = ib.ibObj(f.lCell().centroid());
        auto rCellIbObj = ib.ibObj(f.rCell().centroid());

        if(!lCellIbObj == !rCellIbObj)
        {
            Scalar g = f.rCell().volume()/(f.lCell().volume() + f.rCell().volume());
            u(f) = g*u(f.lCell()) + (1. - g)*u(f.rCell());
        }
        else
            u(f) = rCellIbObj ? ibInterpolate(f.lCell(), f.rCell(), *rCellIbObj): ibInterpolate(f.rCell(), f.lCell(), *lCellIbObj);
    }

    u.setBoundaryFaces();
}
