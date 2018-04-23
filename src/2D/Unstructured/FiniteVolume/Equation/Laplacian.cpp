#include "Laplacian.h"

//namespace fv
//{
//    template<>
//    FiniteVolumeEquation<Vector2D> laplacian(Scalar gamma, FiniteVolumeField<Vector2D> &phi, const CellGroup &cells, Scalar theta)
//    {
//        FiniteVolumeEquation<Vector2D> eqn(phi);
//        const VectorFiniteVolumeField& phi0 = phi.oldField(0);
//
//        for (const Cell &cell: cells)
//        {
//            for (const InteriorLink &nb: cell.neighbours())
//            {
//                Scalar coeff = theta * gamma * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
//                Scalar coeff0 = (1. - theta) * gamma * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
//
//                eqn.add(cell, nb.cell(), coeff);
//                eqn.add(cell, cell, -coeff);
//                eqn.addSource(cell, coeff0 * (phi0(nb.cell()) - phi(cell)));
//
////                Tensor2D tau = theta * gamma * outer(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
////                Tensor2D tau0 = (1. - theta) * gamma * outer(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
////
////                eqn.add(cell, nb.cell(), tau);
////                eqn.add(cell, cell, -tau);
////                eqn.addSource(cell, dot(tau0, phi0(nb.cell()) - phi0(cell)));
//            }
//
//            for (const BoundaryLink &bd: cell.boundaries())
//            {
//                Scalar coeff = theta * gamma * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();
//                Scalar coeff0 = (1. - theta) * gamma * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();
//                Tensor2D tau = theta * gamma * outer(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();
//                Tensor2D tau0 = (1. - theta) * gamma * outer(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();
//
//                switch (phi.boundaryType(bd.face()))
//                {
//                    case VectorFiniteVolumeField::FIXED:
//                        eqn.addSource(cell, coeff * phi(bd.face()));
//                        eqn.add(cell, cell, -coeff);
//                        eqn.addSource(cell, coeff0 * (phi0(bd.face()) - phi0(cell)));
//
////                        eqn.addSource(cell, dot(tau, phi(bd.face())));
////                        eqn.add(cell, cell, -tau);
////                        eqn.addSource(cell, dot(tau0, phi0(bd.face()) - phi0(cell)));
//                        break;
//
//                    case VectorFiniteVolumeField::NORMAL_GRADIENT:
//                    case VectorFiniteVolumeField::SYMMETRY:
//                        break;
//
//                    default:
//                        throw Exception("fv", "VectorFiniteVolumeField", "unrecognized or unspecified boundary type.");
//                }
//            }
//        }
//
//        return eqn;
//    }
//}