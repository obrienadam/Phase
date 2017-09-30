#include "LeeYou.h"
#include "GhostCellImmersedBoundaryObject.h"
#include "CutCell.h"

ScalarFiniteVolumeField lee::massSrc(const VectorFiniteVolumeField &u, const ImmersedBoundary &ib)
{
    ScalarFiniteVolumeField q(u.gridPtr(), "q", 0., false, false);

    for (auto ibObj: ib)
    {
        if (!std::dynamic_pointer_cast<GhostCellImmersedBoundaryObject>(ibObj))
            throw Exception("lee",
                            "massSrc",
                            "mass source can only be computed with \"ghost-cell\" immersed boundary objects.");

        for (const Cell &ibCell: ibObj->ibCells())
            for (const InteriorLink &nb: ibCell.neighbours())
                if (!ibObj->isInIb(nb.cell()))
                {
                    CutCell cell = CutCell(nb.cell(), *ibObj);
                    //q(cell.cell()) -= dot(ibObj->velocity(cell.bFace().center()), cell.solidFaceNorm());

//                    Scalar divS = dot(ibObj->velocity(cell.solidCentroid()), -cell.solidFaceNorm());
//                    for(const CutCellLink& nb: cell.neighbours())
//                    {
//                        if(nb.solidNorm().magSqr() == 0.)
//                            continue;
//
//                        Vector2D us = ibObj->velocity(nb.solidFaceCentroid());
//
//                        const Node& lNode = nb.face().lNode();
//                        const Node& rNode = nb.face().rNode();
//
//                        const Face& f1 = u.grid().interiorFaces().nearestItem(2*lNode - nb.face().centroid());
//                        const Face& f2 = u.grid().interiorFaces().nearestItem(2*rNode - nb.face().centroid());
//
//                        if(ibObj->isInIb(f1) && ibObj->isInIb(f2))
//                        {
//                            divS += dot(us, nb.solidNorm());
//                            continue;
//                        }
//
//                        const Face &f = ibObj->isInIb(f1) ? f2 : f1;
//
//                        Vector2D uf = u(f);
//                        Vector2D uib = ibObj->velocity(nb.xc());
//                        Vector2D eta = (nb.solidFaceCentroid() - f.centroid()).unitVec();
//
//                        Scalar etaF = dot(f.centroid(), eta);
//                        Scalar etaIb = dot(nb.xc(), eta);
//                        Scalar etaS = dot(nb.solidFaceCentroid(), eta);
//
//                        if(etaIb == 0.)
//                        {
//                            uib = ibObj->velocity((f.centroid() + nb.face().centroid()) / 2.);
//                            etaIb = dot((f.centroid() + nb.face().centroid()) / 2., eta);
//                        }
//
//                        Vector2D m = (uib - uf) / (etaIb - etaF);
//                        us = m * (etaS - etaF) + uf;
//                        divS += dot(us, nb.solidNorm());
//                    }
//
//                    q(cell.cell()) = divS;
                }
    }

    return q;
}