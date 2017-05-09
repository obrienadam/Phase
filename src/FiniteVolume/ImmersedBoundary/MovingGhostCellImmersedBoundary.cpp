#include "MovingGhostCellImmersedBoundary.h"
#include "TimeDerivative.h"
#include "FiniteVolumeEquation.h"
#include "CutCell.h"
#include "ForcingCellStencil.h"

Equation<Vector2D> ib::momentumEqn(const ImmersedBoundary& ib,
                                   const ScalarFiniteVolumeField &rho, const ScalarFiniteVolumeField &mu,
                                   ScalarFiniteVolumeField &p, VectorFiniteVolumeField &u, Scalar timeStep)
{
    Equation<Vector2D> uEqn(u, "uEqn");
//    const ScalarFiniteVolumeField &rho0 = rho.prev(0);
//    const VectorFiniteVolumeField &u0 = u.prev(0);

//    for (const Cell &cell: u.grid.cellZone("fluid"))
//    {
//        uEqn.add(cell, cell, cell.volume() * rho(cell) / timeStep);
//        uEqn.addSource(cell, cell.volume() * rho0(cell) * u0(cell) / timeStep);

//        Scalar ac = 0.;
//        for (const InteriorLink &nb: cell.neighbours())
//        {
//            const Face &f = nb.face();
//            const Vector2D &sf = nb.outwardNorm();
//            const Vector2D &rc = nb.rCellVec();

//            Scalar rhoU = rho(f) * dot(u(f), sf);
//            Scalar diff = mu(f) * dot(rc, sf) / dot(rc, rc);

//            Scalar an = std::min(rhoU, 0.) - diff;
//            ac += std::max(rhoU, 0.) + diff;

//            uEqn.add(cell, nb.cell(), an);
//        }

//        for (const BoundaryLink &bd: cell.boundaries())
//        {
//            const Face &f = bd.face();
//            const Vector2D &sf = bd.outwardNorm();
//            const Vector2D &rf = bd.rFaceVec();

//            Scalar rhoU = rho(f)*dot(u(f), sf);
//            Scalar diff = mu(f) * dot(rf, sf) / dot(rf, rf);

//            switch (u.boundaryType(bd.face()))
//            {
//            case VectorFiniteVolumeField::FIXED:
//                uEqn.addSource(cell, (rhoU - diff) * u(f));
//                ac += diff;
//                break;

//            case VectorFiniteVolumeField::NORMAL_GRADIENT:
//                ac += rhoU;
//                break;

//            case VectorFiniteVolumeField::SYMMETRY:
//                break;

//            default:
//                throw Exception("ib", "momentumEqn", "unrecognized or unspecified boundary type.");
//            }
//        }

//        uEqn.add(cell, cell, ac);
//    }

//    for (const ImmersedBoundaryObject &ibObj: ib.ibObjs())
//    {
//        for (const Cell &cell: ibObj.freshlyClearedCells())
//        {
//            uEqn.remove(cell); // Remove equation first
//            ForcingCellStencil st(cell, ibObj.shape(), u.grid.cellZone("fluid"));

//            uEqn.add(st.cell(), st.cell(), 1.);

//            auto cells = st.iCells();
//            auto coeffs = st.iCoeffs();

//            uEqn.add(st.cell(), cells[0], -coeffs[0]);
//            uEqn.add(st.cell(), cells[1], -coeffs[1]);
//            uEqn.addSource(cell, -st.bCoeff()*ibObj.velocity(st.xc()));
//        }

//        for (const GhostCellStencil &stencil: ibObj.ghostCellStencils())
//        {
//            Scalar centralCoeff;
//            std::vector<Scalar> coeffs = stencil.ipCoeffs();

//            //- Boundary assembly
//            switch (ibObj.boundaryType(u.name()))
//            {
//            case ImmersedBoundaryObject::FIXED:
//                centralCoeff = 0.5;
//                for (Scalar &coeff: coeffs)
//                    coeff *= 0.5;

//                uEqn.addSource(stencil.cell(), -ibObj.velocity(stencil.boundaryPoint()));
//                break;

//            default:
//                throw Exception("ib", "momentumEqn", "invalid IB boundary type.");
//            }

//            uEqn.add(stencil.cell(), stencil.cell(), centralCoeff);

//            int i = 0;
//            for (const Cell &ipCell: stencil.ipCells())
//                uEqn.add(stencil.cell(), ipCell, coeffs[i++]);
//        }
//    }

    return uEqn;
}

Equation<Scalar> ib::pressureEqn(const ImmersedBoundary& ib,
                                 const ScalarFiniteVolumeField &rho, const VectorFiniteVolumeField &u,
                                 ScalarFiniteVolumeField &p, Scalar timeStep)
{
    Equation<Scalar> pEqn(p, "pEqn");

//    for (const CutCell& cell: ib.constructCutCells(p.grid.cellZone("fluid")))
//    {
//        Scalar ac = 0.;
//        Scalar divU = 0.;

//        for (const CutCellLink &nb: cell.neighbours())
//        {
//            const Face &f = nb.face();
//            Vector2D sf = nb.fluidNorm();
//            const Vector2D &rc = nb.rCellVec();

//            Scalar diff = timeStep / rho(f) * dot(rc, sf) / dot(rc, rc);

//            Scalar an = diff;
//            ac -= diff;

//            pEqn.add(cell.cell(), nb.cell(), an);
//            divU += dot(u(f), sf);
//        }

//        for (const BoundaryLink &bd: cell.cell().boundaries())
//        {
//            const Face &f = bd.face();
//            const Vector2D &sf = bd.outwardNorm();
//            const Vector2D &rf = bd.rFaceVec();

//            Scalar diff = timeStep / rho(f) * dot(rf, sf) / dot(rf, rf);

//            switch (p.boundaryType(bd.face()))
//            {
//            case ScalarFiniteVolumeField::FIXED:
//                pEqn.addSource(cell.cell(), diff * p(f));
//                ac -= diff;
//                break;

//            case VectorFiniteVolumeField::NORMAL_GRADIENT:
//            case VectorFiniteVolumeField::SYMMETRY:
//                break;

//            default:
//                throw Exception("ib", "pressureEqn", "unrecognized or unspecified boundary type.");
//            }

//            divU += dot(u(f), sf);
//        }

//        if(cell.intersectsIbObj())
//            divU += dot(cell.ibObj().velocity(cell.bFace().center()), -cell.solidFaceNorm());
//        else
//        {
//            if(cell.solidVolume() != 0)
//                throw Exception("ib", "pressureEqn", "bad.");
//        }

//        pEqn.add(cell.cell(), cell.cell(), ac);
//        pEqn.addSource(cell.cell(), -divU);
//    }

//    for (const ImmersedBoundaryObject &ibObj: ib.ibObjs())
//    {
//        for (const GhostCellStencil &stencil: ibObj.ghostCellStencils())
//        {
//            CutCell cutCell(stencil.cell(), ibObj);
//            Scalar divU = 0, sumA = 0;

//            for(const CutCellLink& nb: cutCell.neighbours())
//                if(p.grid.cellZone("fluid").isInGroup(nb.cell()))
//                {
//                    Vector2D n = cutCell.solidFaceNorm().unitVec();
//                    n = Vector2D(n.x*n.x, n.y*n.y);

//                    divU += dot(u(nb.face()), nb.fluidNorm());
//                    sumA += fabs(dot(nb.fluidNorm(), n));
//                }

//            if(cutCell.intersectsIbObj())
//                divU += dot(cutCell.ibObj().velocity(cutCell.bFace().center()), -cutCell.solidFaceNorm());

//            if(sumA > 0) // If sumA is 0, then the IB cell is completely inside the solid
//            {
//                for(const CutCellLink& nb: cutCell.neighbours())
//                    if(p.grid.cellZone("fluid").isInGroup(nb.cell()))
//                    {
//                        Vector2D n = cutCell.solidFaceNorm().unitVec();
//                        n = Vector2D(n.x*n.x, n.y*n.y);

//                        Scalar beta = fabs(dot(nb.fluidNorm(), n))/sumA;
//                        pEqn.addSource(nb.cell(), -beta*divU);
//                    }
//            }
//            else if(cutCell.fluidVolume() > 0)
//            {
//                throw Exception("ibPressureEquation", "pressureEqn", "Ghost-cell has non-zero fluid volume but did not transfer any mass!");
//            }
//            else
//            {
//                if(cutCell.fluidVolume() > 0)
//                {
//                    std::cout << "Cell " << cutCell.cell().id() << " has a non-zero fluid volume, but did not transfer any mass source.\n"
//                              << "alpha = " << cutCell.fluidVolume()/cutCell.totalVolume() << "\n"
//                              << "Solid norm = " << cutCell.solidFaceNorm() << "\n"
//                              << "sumA = " << sumA << std::endl;

//                    for(const CutCellLink& nb: cutCell.neighbours())
//                    {
//                        std::cout << "Fluid norm = " << nb.fluidNorm() << "\n";
//                    }
//                }
//            }

//            Scalar centralCoeff;
//            std::vector<Scalar> coeffs = stencil.ipCoeffs();

//            //- Boundary assembly
//            switch (ibObj.boundaryType(p.name()))
//            {
//            case ImmersedBoundaryObject::NORMAL_GRADIENT:
//                centralCoeff = -1.;
//                break;

//            default:
//                throw Exception("ib", "pressureEqn", "invalid IB boundary type.");
//            }

//            pEqn.add(stencil.cell(), stencil.cell(), centralCoeff);

//            int i = 0;
//            for (const Cell &ipCell: stencil.ipCells())
//                pEqn.add(stencil.cell(), ipCell, coeffs[i++]);
//        }
//    }

    return pEqn;
}

void ib::correctVelocity(const ImmersedBoundary &ib, const ScalarFiniteVolumeField &rho, ScalarFiniteVolumeField &p, VectorFiniteVolumeField& gradP, VectorFiniteVolumeField &u, Scalar timeStep)
{
//    p.setBoundaryFaces();

//    for (const Face &face: u.grid.interiorFaces())
//    {
//        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
//        gradP(face) = (p(face.rCell()) - p(face.lCell()))*rc/dot(rc, rc);
//        u(face) -= timeStep / rho(face) * gradP(face);
//        if(face.lCell().id() == 9763)
//            std::cout << "Cell " << face.lCell().id() << " face vel = " << u(face) << std::endl;
//    }

//    for (const Face &face: u.grid.boundaryFaces())
//    {
//        Vector2D rf = face.centroid() - face.lCell().centroid();
//        gradP(face) = (p(face) - p(face.lCell()))*rf/dot(rf, rf);

//        switch (u.boundaryType(face))
//        {
//        case VectorFiniteVolumeField::FIXED:
//            break;

//        case VectorFiniteVolumeField::SYMMETRY:
//        {
//            const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
//            u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
//        }
//            break;

//        case VectorFiniteVolumeField::NORMAL_GRADIENT:
//            u(face) -= timeStep / rho(face) * gradP(face);
//            break;
//        };
//    }

//    //- Look at trying this with gauss instead!
//    for (const CutCell& cell: ib.constructCutCells(p.grid.cellZone("fluid")))
//    {
//        Scalar sumSfx = 0., sumSfy = 0.;
//        Vector2D &gradPc = gradP(cell.cell()) = Vector2D(0., 0.);

//        for (const CutCellLink &nb: cell.neighbours())
//        {
//            Vector2D sf = nb.fluidNorm();
//            const Vector2D& gradPf = gradP(nb.face());

//            gradPc += Vector2D(gradPf.x * fabs(sf.x), gradPf.y * fabs(sf.y));

//            sumSfx += fabs(sf.x);
//            sumSfy += fabs(sf.y);
//        }

//        for (const BoundaryLink &bd: cell.cell().boundaries())
//        {
//            const Vector2D &sf = bd.outwardNorm();
//            const Vector2D& gradPf = gradP(bd.face());

//            gradPc += Vector2D(gradPf.x * fabs(sf.x), gradPf.y * fabs(sf.y));
//            sumSfx += fabs(sf.x);
//            sumSfy += fabs(sf.y);
//        }

//        const Vector2D &sfs = cell.solidFaceNorm(); // Assuming here that dp/dn = (0,0)

//        sumSfx += fabs(sfs.x);
//        sumSfy += fabs(sfs.y);

//        gradPc = Vector2D(gradPc.x / sumSfx, gradPc.y / sumSfy);
//        u(cell.cell()) -= timeStep/rho(cell.cell())*gradPc;
//    }
}
