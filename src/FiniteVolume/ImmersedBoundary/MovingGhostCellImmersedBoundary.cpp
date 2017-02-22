#include "MovingGhostCellImmersedBoundary.h"

Equation<Vector2D> ib::mv_gc(const std::vector<Ref<const ImmersedBoundaryObject> > &ibObjs, const ScalarFiniteVolumeField& rho, VectorFiniteVolumeField &u, Scalar timeStep)
{
    Equation<Vector2D> uEqn(u);

    for(const ImmersedBoundaryObject& ibObj: ibObjs)
    {
        for(const GhostCellStencil& stencil: ibObj.stencils())
        {
            Scalar centralCoeff;
            std::vector<Scalar> coeffs = stencil.ipCoeffs();

            //- Boundary assembly
            switch(ibObj.boundaryType(u.name()))
            {
            case ImmersedBoundaryObject::FIXED:
                centralCoeff = 0.5;
                for(Scalar &coeff: coeffs)
                    coeff *= 0.5;

                uEqn.addBoundary(stencil.cell(), ibObj.velocity(stencil.boundaryPoint()));

                break;

            default:
                throw Exception("ib", "mv_gc", "invalid boundary type.");
            }

            uEqn.add(stencil.cell(), stencil.cell(), centralCoeff);

            int i = 0;
            for(const Cell& ipCell: stencil.ipCells())
                uEqn.add(stencil.cell(), ipCell, coeffs[i++]);
        }

        // Correct the time integration for freshly cleared cells
        const VectorFiniteVolumeField& u0 = u.prev(0);
        const ScalarFiniteVolumeField& rho0 = rho.prev(0);

        //- This is not perfect, only really works well for cylinders
        for(const Cell& cell: ibObj.freshlyClearedCells())
        {
            Line2D traj(cell.centroid(), ibObj.trajectory().normalVec());
            std::vector<Point2D> intersects = ibObj.shape().intersections(traj);

            Point2D xi = intersects[0];
            for(const Point2D& pt: intersects)
                if((pt - cell.centroid()).magSqr() < (xi - cell.centroid()).magSqr())
                    xi = pt;

            Vector2D uAlpha = ibObj.velocity(xi);
            Scalar timeStepAlpha = (xi - cell.centroid()).mag()/ibObj.trajectory().mag();

            uEqn.add(cell, cell, -rho(cell)*cell.volume()/timeStep);
            uEqn.addBoundary(cell, -rho0(cell)*cell.volume()*u0(cell)/timeStep);

            uEqn.add(cell, cell, rho(cell)*cell.volume()/timeStepAlpha);
            uEqn.addBoundary(cell, rho0(cell)*cell.volume()*uAlpha/timeStepAlpha);
        }

        for(const Cell& cell: ibObj.solidCells())
            u(cell) = ibObj.velocity(cell.centroid());
    }


    return uEqn;
}

Equation<Scalar> ib::mv_gc(const std::vector<Ref<const ImmersedBoundaryObject> > &ibObjs, const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &dp)
{
    Equation<Scalar> pEqn(dp);

    for(const ImmersedBoundaryObject& ibObj: ibObjs)
    {
        for(const GhostCellStencil& stencil: ibObj.stencils())
        {
            Scalar centralCoeff, l = stencil.length();
            std::vector<Scalar> coeffs = stencil.ipCoeffs();

            //- Boundary assembly
            switch(ibObj.boundaryType(dp.name()))
            {

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                centralCoeff = -1.;
                break;

            default:
                throw Exception("ib", "mv_gc", "invalid boundary type.");
            }

            pEqn.add(stencil.cell(), stencil.cell(), centralCoeff);

            int i = 0;
            for(const Cell& ipCell: stencil.ipCells())
                pEqn.add(stencil.cell(), ipCell, coeffs[i++]);
        }
    }

    return pEqn;
}
