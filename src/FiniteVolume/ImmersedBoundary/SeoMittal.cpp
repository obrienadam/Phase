#include "SeoMittal.h"
#include "GhostCellImmersedBoundaryObject.h"

Equation<Scalar> seo::laplacian(const ImmersedBoundary &ib, Scalar rho, Scalar timeStep, ScalarFiniteVolumeField& p)
{
    Equation<Scalar> pEqn(p, "pEqn");

    for(const CutCell &cell: ib.constructCutCells(ib.zone()))
    {
        Scalar ac = 0.;

        for(const CutCellLink &nb: cell.neighbours())
        {
            const Vector2D& rc = nb.rCellVec();
            Vector2D sf = nb.fluidNorm();

            Scalar anb = timeStep/rho*dot(rc, sf)/dot(rc, rc);
            ac -= anb;

            pEqn.add(cell, nb.cell(), anb);
        }

        for(const BoundaryLink& bd: cell.boundaries())
        {
            const Face &face = bd.face();
            const Vector2D &sf = bd.outwardNorm();
            const Vector2D &rf = bd.rFaceVec();
            Scalar anb = timeStep/rho*dot(rf, sf)/dot(rf, rf);

            switch(p.boundaryType(face))
            {
                case ScalarFiniteVolumeField::FIXED:
                    pEqn.addSource(cell, anb*p(face));
                    ac -= anb;
                    break;
                case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                case ScalarFiniteVolumeField::SYMMETRY:
                    break;
                default:
                    throw Exception("seo", "laplacian", "invalid boundary type.");
            }
        }

        pEqn.add(cell, cell, ac);

        if(cell.intersectsIbObj())
        {
            Vector2D sf = cell.solidFaceNorm();
            Vector2D dP = -rho*cell.ibObj().acceleration(cell.bFace().center());
            Scalar flux = timeStep/rho*dot(dP, sf);
            pEqn.addSource(cell, flux);
        }
    }

    for(const ImmersedBoundaryObject& ibObj: ib.ibObjs())
    {
        try
        {
            const auto &gcIbObj = dynamic_cast<const GhostCellImmersedBoundaryObject &>(ibObj);

            for (const GhostCellStencil &st: gcIbObj.stencils())
            {
                Scalar dPdN = dot(-rho * ibObj.acceleration(st.boundaryPoint()), st.unitNormal());

                pEqn.add(st.cell(), st.cell(), 1. / st.length());
                pEqn.addSource(st.cell(), -dPdN);

                std::vector<Scalar> coeffs = st.ipCoeffs();
                int i = 0;
                for (const Cell &ipCell: st.ipCells())
                    pEqn.add(st.cell(), ipCell, -coeffs[i++] / st.length());
            }
        }
        catch (std::bad_cast e)
        {
            throw Exception("seo", "laplacian", "failed to cast ImmersedBoundaryObject to GhostCellImmersedBoundaryObject.");
        }
    }

    return pEqn;
}

ScalarFiniteVolumeField seo::div(const ImmersedBoundary &ib, const VectorFiniteVolumeField &u)
{
    ScalarFiniteVolumeField divU(u.grid(), "divU", 0);

    for (const CutCell &cell: ib.constructCutCells(ib.zone()))
    {
        for (const CutCellLink &nb: cell.neighbours())
            divU(cell) += dot(u(nb.face()), nb.fluidNorm());

        for (const BoundaryLink &bd: cell.boundaries())
            divU(cell) += dot(u(bd.face()), bd.outwardNorm());

        if(cell.intersectsIbObj())
            divU(cell) += dot(cell.ibObj().velocity(cell.bFace().center()), cell.solidFaceNorm());
    }

    for (const ImmersedBoundaryObject &ibObj: ib.ibObjs())
        for (const CutCell &cell: ib.constructCutCells(ibObj.ibCells()))
        {
            for (const CutCellLink &nb: cell.neighbours())
                divU(cell) += dot(u(nb.face()), nb.fluidNorm());

            for (const BoundaryLink &bd: cell.boundaries())
                divU(cell) += dot(u(bd.face()), bd.outwardNorm());

            if(cell.intersectsIbObj())
                divU(cell) += dot(ibObj.velocity(cell.bFace().center()), cell.solidFaceNorm());

            Scalar sumBeta = 0.;
            std::vector<Scalar> beta;

            for (const CutCellLink &nb: cell.neighbours())
                if (!ibObj.isInIb(nb.cell().centroid()))
                {
                    beta.push_back(fabs(dot(nb.fluidNorm(), cell.solidFaceNorm())));
                    sumBeta += beta.back();
                }

            int i = 0;
            if(sumBeta > 0)
            {
                for (const CutCellLink &nb: cell.neighbours())
                    if (!ibObj.isInIb(nb.cell().centroid()))
                        divU(nb.cell()) += beta[i++] / sumBeta * divU(cell);
            }

            divU(cell) = 0.;
        }

    return divU;
}

void seo::correct(const ImmersedBoundary &ib, Scalar rho, const ScalarFiniteVolumeField &p, VectorFiniteVolumeField &gradP,
                  VectorFiniteVolumeField &u, Scalar timeStep)
{
    //- Correct faces
    for (const Face &face: u.grid().interiorFaces())
    {
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        gradP(face) = (p(face.rCell()) - p(face.lCell())) * rc / dot(rc, rc);
        u(face) -= timeStep / rho * gradP(face);
    }

    for(const Face& face: u.grid().boundaryFaces())
    {
        Vector2D rf = face.centroid() - face.lCell().centroid();
        gradP(face) = (p(face) - p(face.lCell())) * rf / dot(rf, rf);
    }

    for(const Patch& patch: u.grid().patches())
    {
        switch(u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                for(const Face& face: patch)
                {
                    Vector2D nWall = face.outwardNorm(face.lCell().centroid());
                    u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
                }
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for(const Face &face: patch)
                    u(face) -= timeStep / rho * gradP(face);
                break;
        }
    }

    //- Correct cell center
    for (const CutCell &cell: ib.constructCutCells(ib.zone()))
    {
        Scalar sumAx = 0., sumAy = 0.;
        gradP(cell) = Vector2D(0., 0.);

        for(const CutCellLink& nb: cell.neighbours())
        {
            Vector2D sf = nb.fluidNorm();
            sf = Vector2D(fabs(sf.x), fabs(sf.y));

            sumAx += sf.x;
            sumAy += sf.y;

            gradP(cell) += Vector2D(gradP(nb.face()).x*sf.x, gradP(nb.face()).y*sf.y);
        }

        for(const BoundaryLink& bd: cell.boundaries())
        {
            Scalar ax = fabs(bd.outwardNorm().x), ay = fabs(bd.outwardNorm().y);

            sumAx += ax;
            sumAy += ay;

            gradP(cell) += Vector2D(gradP(bd.face()).x*ax, gradP(bd.face()).y*ay);
        }

        if(cell.intersectsIbObj())
        {
            Scalar ax = fabs(cell.solidFaceNorm().x), ay = fabs(cell.solidFaceNorm().y);

            sumAx += ax;
            sumAy += ay;

            Vector2D sf = cell.solidFaceNorm();
            Vector2D dP = -rho*cell.ibObj().acceleration(cell.bFace().center());
            dP = dot(dP, sf)*sf/sf.magSqr();

            gradP(cell) += Vector2D(dP.x*ax, dP.y*ay);
        }

        gradP(cell) = Vector2D(gradP(cell).x/sumAx, gradP(cell).y/sumAy);
        u(cell) -= timeStep/rho*gradP(cell);
    }
}

Equation<Scalar> seo::laplacian(const ImmersedBoundary &ib, const ScalarFiniteVolumeField& rho, Scalar timeStep,
                                ScalarFiniteVolumeField &p)
{
    Equation<Scalar> pEqn(p, "pEqn");

    for(const CutCell &cell: ib.constructCutCells(ib.zone()))
    {
        Scalar ac = 0.;

        for(const CutCellLink &nb: cell.neighbours())
        {
            const Vector2D& rc = nb.rCellVec();
            Vector2D sf = nb.fluidNorm();

            Scalar anb = timeStep/rho(nb.face())*dot(rc, sf)/dot(rc, rc);
            ac -= anb;

            pEqn.add(cell, nb.cell(), anb);
        }

        for(const BoundaryLink& bd: cell.boundaries())
        {
            const Face &face = bd.face();
            const Vector2D &sf = bd.outwardNorm();
            const Vector2D &rf = bd.rFaceVec();
            Scalar anb = timeStep/rho(bd.face())*dot(rf, sf)/dot(rf, rf);

            switch(p.boundaryType(face))
            {
                case ScalarFiniteVolumeField::FIXED:
                    pEqn.addSource(cell, anb*p(face));
                    ac -= anb;
                    break;
                case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                case ScalarFiniteVolumeField::SYMMETRY:
                    break;
                default:
                    throw Exception("seo", "laplacian", "invalid boundary type.");
            }
        }

        pEqn.add(cell, cell, ac);

        if(cell.intersectsIbObj())
        {
            Vector2D sf = cell.solidFaceNorm();
            Vector2D dP = -rho(cell)*cell.ibObj().acceleration(cell.bFace().center());
            Scalar flux = timeStep/rho(cell)*dot(dP, sf);
            pEqn.addSource(cell, flux);
        }
    }

    for(const ImmersedBoundaryObject& ibObj: ib.ibObjs())
    {
        try
        {
            const auto &gcIbObj = dynamic_cast<const GhostCellImmersedBoundaryObject &>(ibObj);

            for (const GhostCellStencil &st: gcIbObj.stencils())
            {
                Scalar rhoF = st.ipValue(rho);
                Scalar dPdN = dot(-rhoF * ibObj.acceleration(st.boundaryPoint()), st.unitNormal());

                pEqn.add(st.cell(), st.cell(), 1. / st.length());
                pEqn.addSource(st.cell(), -dPdN);

                std::vector<Scalar> coeffs = st.ipCoeffs();
                int i = 0;
                for (const Cell &ipCell: st.ipCells())
                    pEqn.add(st.cell(), ipCell, -coeffs[i++] / st.length());
            }
        }
        catch (std::bad_cast e)
        {
            throw Exception("seo", "laplacian", "failed to cast ImmersedBoundaryObject to GhostCellImmersedBoundaryObject.");
        }
    }

    return pEqn;
}


void seo::correct(const ImmersedBoundary &ib, const ScalarFiniteVolumeField &rho, const ScalarFiniteVolumeField& p,
                  VectorFiniteVolumeField &gradP, VectorFiniteVolumeField &u, Scalar timeStep)
{
    //- Correct faces
    for (const Face &face: u.grid().interiorFaces())
    {
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        gradP(face) = (p(face.rCell()) - p(face.lCell())) * rc / dot(rc, rc);
        u(face) -= timeStep / rho(face) * gradP(face);
    }

    for(const Face& face: u.grid().boundaryFaces())
    {
        Vector2D rf = face.centroid() - face.lCell().centroid();
        gradP(face) = (p(face) - p(face.lCell())) * rf / dot(rf, rf);
    }

    for(const Patch& patch: u.grid().patches())
    {
        switch(u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                for(const Face& face: patch)
                {
                    Vector2D nWall = face.outwardNorm(face.lCell().centroid());
                    u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
                }
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for(const Face &face: patch)
                    u(face) -= timeStep / rho(face) * gradP(face);
                break;
        }
    }

    //- Correct cell center
    for (const CutCell &cell: ib.constructCutCells(ib.zone()))
    {
        Scalar sumAx = 0., sumAy = 0.;
        gradP(cell) = Vector2D(0., 0.);

        for(const CutCellLink& nb: cell.neighbours())
        {
            Vector2D sf = nb.fluidNorm();
            sf = Vector2D(fabs(sf.x), fabs(sf.y));

            sumAx += sf.x;
            sumAy += sf.y;

            gradP(cell) += Vector2D(gradP(nb.face()).x*sf.x, gradP(nb.face()).y*sf.y);
        }

        for(const BoundaryLink& bd: cell.boundaries())
        {
            Scalar ax = fabs(bd.outwardNorm().x), ay = fabs(bd.outwardNorm().y);

            sumAx += ax;
            sumAy += ay;

            gradP(cell) += Vector2D(gradP(bd.face()).x*ax, gradP(bd.face()).y*ay);
        }

        if(cell.intersectsIbObj())
        {
            Scalar ax = fabs(cell.solidFaceNorm().x), ay = fabs(cell.solidFaceNorm().y);

            sumAx += ax;
            sumAy += ay;

            Vector2D sf = cell.solidFaceNorm();
            Vector2D dP = -rho(cell)*cell.ibObj().acceleration(cell.bFace().center());
            dP = dot(dP, sf)*sf/sf.magSqr();

            gradP(cell) += Vector2D(dP.x*ax, dP.y*ay);
        }

        gradP(cell) = Vector2D(gradP(cell).x/sumAx, gradP(cell).y/sumAy);
        u(cell) -= timeStep/rho(cell)*gradP(cell);
    }
}

Scalar seo::maxDivergence(const ImmersedBoundary &ib, const VectorFiniteVolumeField& u)
{
    Scalar maxDivU = 0.;

    for (const CutCell &cell: ib.constructCutCells(ib.zone()))
    {
        Scalar divU = 0.;

        for (const CutCellLink &nb: cell.neighbours())
            divU += dot(u(nb.face()), nb.fluidNorm());

        for (const BoundaryLink &bd: cell.boundaries())
            divU += dot(u(bd.face()), bd.outwardNorm());

        if(cell.intersectsIbObj())
            divU += dot(cell.ibObj().velocity(cell.bFace().center()), cell.solidFaceNorm());

        if(fabs(divU) > fabs(maxDivU))
            maxDivU = divU;
    }

    return maxDivU;
}