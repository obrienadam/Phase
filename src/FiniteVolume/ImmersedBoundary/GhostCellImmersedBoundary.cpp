#include "GhostCellImmersedBoundary.h"
#include "BilinearInterpolation.h"

namespace gc {

Equation<ScalarFiniteVolumeField> ib(const std::vector<ImmersedBoundaryObject>& ibObjs, ScalarFiniteVolumeField& field)
{
    Equation<ScalarFiniteVolumeField> eqn(field);

    for(const ImmersedBoundaryObject& ibObj: ibObjs)
        for(const Cell &cell: ibObj.cells())
        {
            const Index row = cell.localIndex();
            const Point2D& imagePoint = ibObj.imagePoint(cell);

            const std::vector< Ref<const Cell> > &kNN = ibObj.imagePointCells(cell);
            const BilinearInterpolation &bi = ibObj.imagePointInterpolation(cell);

            std::vector<Scalar> coeffs = bi(imagePoint);

            std::vector<Index> cols = {
                kNN[0].get().localIndex(),
                kNN[1].get().localIndex(),
                kNN[2].get().localIndex(),
                kNN[3].get().localIndex(),
            };

            Scalar centralCoeff, l;
            // Then here we shall do the boundary assembly!
            switch(ibObj.boundaryType(field.name()))
            {
            case ImmersedBoundaryObject::FIXED:
                centralCoeff = 1.;
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                centralCoeff = -1.;
                break;

            case ImmersedBoundaryObject::CONTACT_ANGLE:
                centralCoeff = -1;
                break;

            case ImmersedBoundaryObject::PARTIAL_SLIP:
                l = (cell.centroid() - imagePoint).mag();
                centralCoeff = 1. + 2.*ibObj.boundaryRefValue(field.name())/l;

                for(Scalar &coeff: coeffs)
                    coeff *= 1. - 2.*ibObj.boundaryRefValue(field.name())/l;

                break;

            default:
                throw Exception("gc", "ib", "invalid boundary type.");
            }

            eqn.add(row, row, centralCoeff);

            for(int i = 0; i < coeffs.size(); ++i)
                eqn.add(row, cols[i], coeffs[i]);
        }

    return eqn;
}

Equation<VectorFiniteVolumeField> ib(const std::vector<ImmersedBoundaryObject> &ibObjs, VectorFiniteVolumeField &field)
{
    Equation<VectorFiniteVolumeField> eqn(field);
    const Size nActiveCells = field.grid.nActiveCells();

    for(const ImmersedBoundaryObject &ibObj: ibObjs)
        for(const Cell &cell: ibObj.cells())
        {
            const Index rowX = cell.localIndex();
            const Index rowY = rowX + nActiveCells;

            const Point2D imagePoint = ibObj.imagePoint(cell);

            const std::vector< Ref<const Cell> > &kNN = ibObj.imagePointCells(cell);
            const BilinearInterpolation &bi = ibObj.imagePointInterpolation(cell);

            std::vector<Scalar> coeffsX = bi(imagePoint), coeffsY;
            coeffsY = coeffsX;

            std::vector<Index> colsX = {
                kNN[0].get().localIndex(),
                kNN[1].get().localIndex(),
                kNN[2].get().localIndex(),
                kNN[3].get().localIndex(),
            };

            std::vector<Index> colsY(4);
            for(int i = 0; i < 4; ++i)
                colsY[i] = colsX[i] + nActiveCells;

            // Then here we shall do the boundary assembly!
            Scalar centralCoeffX, centralCoeffY, l, tmpScalar;
            Vector2D n, t;
            switch(ibObj.boundaryType(field.name()))
            {
            case ImmersedBoundaryObject::FIXED:
                centralCoeffX = centralCoeffY = 1.;
                break;
            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                centralCoeffX = centralCoeffY = -1.;
                break;

            case ImmersedBoundaryObject::CONTACT_ANGLE:
                throw Exception("gc", "ib", "contact angle boundaries are not valid for vector fields.");

            case ImmersedBoundaryObject::PARTIAL_SLIP:

                //- Note: the tangent equation replaces the u.x equation and the normal equation replaces the u.y equation

                //- Enforce the tangent condition
                n = (imagePoint - cell.centroid()).unitVec();
                t = n.tangentVec();

                l = (imagePoint - cell.centroid()).mag();

                tmpScalar = ibObj.boundaryRefValue(field.name());
                tmpScalar = (1. + 2.*tmpScalar/l)/(1. - 2.*tmpScalar/l);

                centralCoeffX = tmpScalar*t.x;

                colsX.push_back(rowY);
                coeffsX.push_back(tmpScalar*t.y);

                for(int i = 0; i < 4; ++i)
                {
                    colsX.push_back(colsY[i]);
                    coeffsX.push_back(coeffsY[i]*t.y);
                }

                //- Enforce the normal condition
                centralCoeffY = n.y;

                colsY.push_back(rowX);
                coeffsY.push_back(n.x);

                for(int i = 0; i < 4; ++i)
                {
                    colsY.push_back(colsX[i]);
                    coeffsY.push_back(coeffsX[i]*n.x);
                }

                //- Modify the original coefficients
                for(int i = 0; i < 4; ++i)
                {
                    coeffsX[i] *= t.x;
                    coeffsY[i] *= n.y;
                }

                break;

            default:
                throw Exception("gc", "ib", "invalid boundary type.");
            }

            eqn.add(rowX, rowX, centralCoeffX);
            eqn.add(rowY, rowY, centralCoeffY);

            for(int i = 0; i < coeffsX.size(); ++i)
            {
                eqn.add(rowX, colsX[i], coeffsX[i]);
                eqn.add(rowY, colsY[i], coeffsY[i]);
            }
        }

    return eqn;
}


Equation<VectorFiniteVolumeField> ib(const std::vector<ImmersedBoundaryObject> &ibObjs, const ScalarFiniteVolumeField &rho, VectorFiniteVolumeField &field)
{
    Equation<VectorFiniteVolumeField> eqn(field);
    const Size nActiveCells = field.grid.nActiveCells();
    Scalar rhoIp;

    for(const ImmersedBoundaryObject &ibObj: ibObjs)
        for(const Cell &cell: ibObj.cells())
        {
            const Index rowX = cell.localIndex();
            const Index rowY = rowX + nActiveCells;

            const Point2D imagePoint = ibObj.imagePoint(cell);

            const std::vector< Ref<const Cell> > &kNN = ibObj.imagePointCells(cell);
            const BilinearInterpolation &bi = ibObj.imagePointInterpolation(cell);

            std::vector<Scalar> coeffsX = bi(imagePoint), coeffsY;
            coeffsY = coeffsX;

            std::vector<Index> colsX = {
                kNN[0].get().localIndex(),
                kNN[1].get().localIndex(),
                kNN[2].get().localIndex(),
                kNN[3].get().localIndex(),
            };

            std::vector<Scalar> rhoVals = {
                rho(kNN[0]),
                rho(kNN[1]),
                rho(kNN[2]),
                rho(kNN[3])
            };

            rhoIp = bi(rhoVals, imagePoint);

            std::vector<Index> colsY(4);
            for(int i = 0; i < 4; ++i)
                colsY[i] = colsX[i] + nActiveCells;

            // Then here we shall do the boundary assembly!
            Scalar centralCoeffX, centralCoeffY, l, tmpScalar;
            Vector2D n, t;
            switch(ibObj.boundaryType(field.name()))
            {
            case ImmersedBoundaryObject::FIXED:
                centralCoeffX = centralCoeffY = 1.;
                break;
            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                centralCoeffX = centralCoeffY = -1.;
                break;

            case ImmersedBoundaryObject::CONTACT_ANGLE:
                throw Exception("gc", "ib", "contact angle boundaries are not valid for vector fields.");

            case ImmersedBoundaryObject::PARTIAL_SLIP:

                //- Note: the tangent equation replaces the u.x equation and the normal equation replaces the u.y equation

                //- Enforce the tangent condition
                n = (imagePoint - cell.centroid()).unitVec();
                t = n.tangentVec();

                l = (imagePoint - cell.centroid()).mag();

                tmpScalar = ibObj.boundaryRefValue(field.name());
                tmpScalar = (1. + 2.*tmpScalar/l)/(1. - 2.*tmpScalar/l);

                centralCoeffX = tmpScalar*t.x;

                colsX.push_back(rowY);
                coeffsX.push_back(tmpScalar*t.y);

                for(int i = 0; i < 4; ++i)
                {
                    colsX.push_back(colsY[i]);
                    coeffsX.push_back(coeffsY[i]*t.y);
                }

                //- Enforce the normal condition
                centralCoeffY = n.y;

                colsY.push_back(rowX);
                coeffsY.push_back(n.x);

                for(int i = 0; i < 4; ++i)
                {
                    colsY.push_back(colsX[i]);
                    coeffsY.push_back(coeffsX[i]*n.x);
                }

                //- Modify the original coefficients
                for(int i = 0; i < 4; ++i)
                {
                    coeffsX[i] *= t.x;
                    coeffsY[i] *= n.y;
                }

                break;

            default:
                throw Exception("gc", "ib", "invalid boundary type.");
            }

            eqn.add(rowX, rowX,  centralCoeffX*rhoIp/(rho(cell) + rhoIp));
            eqn.add(rowY, rowY, centralCoeffY*rhoIp/(rho(cell) + rhoIp));

            for(int i = 0; i < coeffsX.size(); ++i)
            {
                eqn.add(rowX, colsX[i], coeffsX[i]*rho(cell)/(rho(cell) + rhoIp));
                eqn.add(rowY, colsY[i], coeffsY[i]*rho(cell)/(rho(cell) + rhoIp));
            }
        }

    return eqn;
}

}
