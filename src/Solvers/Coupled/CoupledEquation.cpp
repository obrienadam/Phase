#include "CoupledEquation.h"

CoupledEquation::CoupledEquation(const ScalarFiniteVolumeField &rho, const ScalarFiniteVolumeField &mu, VectorFiniteVolumeField &u, ScalarFiniteVolumeField &p)
    :
      rho_(rho),
      mu_(mu),
      u_(u),
      gradP_(p.grid, "gradP"),
      p_(p),
      d_(u_.grid, "d")
{
    nActiveCells_ = u.grid.nActiveCells();
    nVars_ = 3*nActiveCells_;

    spMat_ = SparseMatrix(nVars_, nVars_, 13);
    rhs_ = SparseVector::Zero(nVars_);
}

Scalar CoupledEquation::solve(Scalar timeStep)
{
    spMat_.setZero();
    rhs_.setZero();

    u_.savePreviousTimeStep(timeStep, 1);

    printf("Assembling momentum equation...\n");
    rhieChowInterpolation();
    assembleMomentumEquation(timeStep);

    printf("Assembling continuity equation...\n");
    computeD();
    assembleContinuityEquation();

    printf("Solving coupled Navier-Stokes equations...\n");
    SparseVector x = spMat_.solve(rhs_);

    for(const Cell &cell: u_.grid.activeCells())
        u_[cell.id()].x = x(cell.globalIndex());

    for(const Cell &cell: u_.grid.activeCells())
        u_[cell.id()].y = x(cell.globalIndex() + nActiveCells_);

    for(const Cell &cell: p_.grid.activeCells())
        p_[cell.id()] = x(cell.globalIndex() + 2*nActiveCells_);

    printf("Solved coupled Navier-Stokes equations. Error = %lf, number of iterations = %d.\n", spMat_.error(), spMat_.nIterations());

    return spMat_.error();
}

//- Protected methods

void CoupledEquation::computeD()
{
    const auto diag = spMat_.diagonal();

    for(const Cell &cell: u_.grid.fluidCells())
        d_[cell.id()] = cell.volume()/diag[cell.globalIndex()];

    interpolateFaces(d_);
}

void CoupledEquation::assembleMomentumEquation(Scalar timeStep)
{
    for(const Cell &cell: u_.grid.fluidCells())
    {
        const size_t rowU = cell.globalIndex();
        const size_t rowV = rowU + nActiveCells_;
        const size_t rowP = rowV + nActiveCells_;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const size_t colU = nb.cell().globalIndex();
            const size_t colV = colU + nActiveCells_;
            const size_t colP = colV + nActiveCells_;

            const Face &face = nb.face();
            const Vector2D& sf = nb.outwardNorm();
            const Vector2D& rc = nb.rCellVec();

            spMat_.coeffRef(rowU, rowU) += cell.volume()/timeStep;
            spMat_.coeffRef(rowV, rowV) += cell.volume()/timeStep;

            rhs_(rowU) += cell.volume()*u_.prev(0)[cell.id()].x/timeStep;
            rhs_(rowV) += cell.volume()*u_.prev(0)[cell.id()].y/timeStep;

            const Scalar massFlux = rho_.faces()[face.id()]*dot(u_.faces()[face.id()], sf);

            spMat_.coeffRef(rowU, rowU) += std::max(massFlux, 0.);
            spMat_.coeffRef(rowV, rowV) += std::max(massFlux, 0.);
            spMat_.coeffRef(rowU, colU) += std::min(massFlux, 0.);
            spMat_.coeffRef(rowV, colV) += std::min(massFlux, 0.);

            const Scalar diffFlux = mu_.faces()[face.id()]*dot(rc, sf)/dot(rc, rc);

            spMat_.coeffRef(rowU, rowU) += diffFlux;
            spMat_.coeffRef(rowV, rowV) += diffFlux;
            spMat_.coeffRef(rowU, colU) -= diffFlux;
            spMat_.coeffRef(rowV, colV) -= diffFlux;

            const Scalar g = nb.cell().volume()/(cell.volume() + nb.cell().volume());

            //- Pressure source
            spMat_.coeffRef(rowU, rowP) += g*sf.x;
            spMat_.coeffRef(rowV, rowP) += g*sf.y;
            spMat_.coeffRef(rowU, colP) += (1. - g)*sf.x;
            spMat_.coeffRef(rowV, colP) += (1. - g)*sf.y;
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Face &face = bd.face();
            const Vector2D& sf = bd.outwardNorm();
            const Vector2D& rf = bd.rFaceVec();

            const Scalar massFlux = rho_.faces()[face.id()]*dot(u_.faces()[face.id()], sf);
            const Scalar diffFlux = mu_.faces()[face.id()]*dot(rf, sf)/dot(rf, rf);

            switch(u_.boundaryType(bd.face().id()))
            {
            case VectorFiniteVolumeField::FIXED:
                spMat_.coeffRef(rowU, rowU) += diffFlux;
                spMat_.coeffRef(rowV, rowV) += diffFlux;

                rhs_(rowU) += (diffFlux - massFlux)*u_.faces()[bd.face().id()].x;
                rhs_(rowV) += (diffFlux - massFlux)*u_.faces()[bd.face().id()].y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                spMat_.coeffRef(rowU, rowU) += massFlux;
                spMat_.coeffRef(rowV, rowV) += massFlux;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                //break;

            default:
                throw Exception("CoupledEquation", "assembleMomentumEquation", "unrecognized or unspecified boundary type.");
            }

            switch(p_.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                rhs_(rowU) -= p_.faces()[bd.face().id()]*sf.x;
                rhs_(rowV) -= p_.faces()[bd.face().id()]*sf.y;
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT: case ScalarFiniteVolumeField::SYMMETRY:
                spMat_.coeffRef(rowU, rowP) += sf.x;
                spMat_.coeffRef(rowV, rowP) += sf.y;
                break;

            default:
                throw Exception("CoupledEquation", "assembleMomentumEquation", "unrecognized or unspecified boundary type.");
            }
        }
    }
}

void CoupledEquation::assembleContinuityEquation()
{

    for(const Cell &cell: p_.grid.fluidCells())
    {
        const size_t rowU = cell.globalIndex();
        const size_t rowV = rowU + nActiveCells_;
        const size_t rowP = rowV + nActiveCells_;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const size_t colU = nb.cell().globalIndex();
            const size_t colV = colU + nActiveCells_;
            const size_t colP = colV + nActiveCells_;

            const Face &face = nb.face();
            const Vector2D& sf = nb.outwardNorm();
            const Vector2D& rc = nb.rCellVec();
            const Scalar df = d_.faces()[face.id()];

            const Scalar diffFlux = df*dot(rc, sf)/dot(rc, rc);

            spMat_.coeffRef(rowP, rowP) -= diffFlux;
            spMat_.coeffRef(rowP, colP) += diffFlux;

            const Scalar g = nb.cell().volume()/(cell.volume() + nb.cell().volume());

            spMat_.coeffRef(rowP, rowU) -= g*sf.x;
            spMat_.coeffRef(rowP, rowV) -= g*sf.y;
            spMat_.coeffRef(rowP, colU) -= (1. - g)*sf.x;
            spMat_.coeffRef(rowP, colV) -= (1. - g)*sf.y;

            rhs_(rowP) += df*dot(gradP_.faces()[face.id()], sf);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D &sf = bd.outwardNorm();
            const Scalar massFlux = dot(u_.faces()[bd.face().id()], sf);
            const Scalar diffFlux = d_.faces()[bd.face().id()]*dot(bd.rFaceVec(), sf)/dot(bd.rFaceVec(), bd.rFaceVec());

            switch(u_.boundaryType(bd.face().id()))
            {
            case VectorFiniteVolumeField::FIXED:
                rhs_(rowP) -= massFlux;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                spMat_.coeffRef(rowP, rowU) -= sf.x;
                spMat_.coeffRef(rowP, rowV) -= sf.y;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                //break;
            default:
                throw Exception("CoupledEquation", "assembleContinuity", "unrecognized or unspecified boundary type.");
            }

            switch(p_.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                rhs_(rowP) -= diffFlux*p_.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:

                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                //break;

            default:
                throw Exception("CoupledEquation", "assembleContinuity", "unrecognized or unspecified boundary type.");
            }
        }
    }
}

void CoupledEquation::rhieChowInterpolation()
{
    interpolateFaces(p_);
    gradP_ = grad(p_);

    interpolateFaces(gradP_);
    interpolateFaces(u_);

    for(const Face& face: u_.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        Vector2D sf = face.outwardNorm(lCell.centroid());
        Vector2D rc = rCell.centroid() - lCell.centroid();

        u_.faces()[face.id()] -= d_.faces()[face.id()]*((p_[rCell.id()] - p_[lCell.id()])*sf/dot(rc, rc) - gradP_.faces()[face.id()]);
    }

    for(const Face& face: u_.grid.boundaryFaces())
    {
        Vector2D nWall;

        switch(u_.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u_.faces()[face.id()] = u_[face.lCell().id()];
            break;

        case VectorFiniteVolumeField::SYMMETRY:
            nWall = face.outwardNorm(face.lCell().centroid()).unitVec();
            u_.faces()[face.id()] = u_[face.lCell().id()] - dot(u_[face.lCell().id()], nWall)*nWall;
            break;

        default:
            throw Exception("CoupledEquation", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}
