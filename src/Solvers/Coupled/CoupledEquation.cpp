#include "CoupledEquation.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

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
    phi_ = SparseVector::Zero(nVars_);

    spMat_.setTolerance(1e-10);
    spMat_.setMaxIterations(1000);
    spMat_.setFill(2);
}

Scalar CoupledEquation::solve(Scalar timeStep)
{
    spMat_.setZero();
    rhs_ = SparseVector::Zero(nVars_);

    u_.savePreviousTimeStep(timeStep, 1);

    printf("Assembling momentum equation...\n");
    rhieChowInterpolation();
    assembleMomentumEquation(timeStep);

    printf("Assembling continuity equation...\n");
    assembleContinuityEquation();

    printf("Assembling global coupled matrix...\n");
    spMat_.assemble(triplets_);
    triplets_.clear();

    printf("Solving coupled Navier-Stokes equations...\n");
    phi_ = spMat_.solve(rhs_, phi_, SparseMatrix::IncompleteLUT, true);
    printf("Solved coupled Navier-Stokes equations. Error = %lf, number of iterations = %d.\n", spMat_.error(), spMat_.nIterations());

    for(const Cell &cell: u_.grid.activeCells())
        u_(cell).x = phi_(cell.globalIndex());

    for(const Cell &cell: u_.grid.activeCells())
        u_(cell).y = phi_(cell.globalIndex() + nActiveCells_);

    for(const Cell &cell: p_.grid.activeCells())
        p_(cell) = phi_(cell.globalIndex() + 2*nActiveCells_);

    return spMat_.error();
}

//- Protected methods

void CoupledEquation::assembleMomentumEquation(Scalar timeStep)
{
    d_.fill(0.);

    for(const Cell &cell: u_.grid.fluidCells())
    {
        const Index rowU = cell.globalIndex();
        const Index rowV = rowU + nActiveCells_;
        const Index rowP = rowV + nActiveCells_;
        Scalar centralCoeffUV = 0.;
        Scalar centralCoeffUP = 0.;
        Scalar centralCoeffVP = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Index colU = nb.cell().globalIndex();
            const Index colV = colU + nActiveCells_;
            const Index colP = colV + nActiveCells_;
            Scalar coeffUV = 0.;
            Scalar coeffUP = 0.;
            Scalar coeffVP = 0.;

            const Face &face = nb.face();
            const Vector2D& sf = nb.outwardNorm();
            const Vector2D& rc = nb.rCellVec();

            //- Time term
            centralCoeffUV += cell.volume()/timeStep;

            rhs_(rowU) += cell.volume()*u_.prev(0)(cell).x/timeStep;
            rhs_(rowV) += cell.volume()*u_.prev(0)(cell).y/timeStep;

            //- Advection term
            const Scalar massFlux = rho_(face)*dot(u_(face), sf);

            centralCoeffUV += std::max(massFlux, 0.);
            coeffUV += std::min(massFlux, 0.);

            //- Diffusion term
            const Scalar diffFlux = mu_(face)*dot(rc, sf)/dot(rc, rc);

            centralCoeffUV += diffFlux;
            coeffUV -= diffFlux;

            //- Pressure term
            const Scalar g = nb.cell().volume()/(cell.volume() + nb.cell().volume());

            centralCoeffUP += g*sf.x;
            centralCoeffVP += g*sf.y;

            coeffUP += (1. - g)*sf.x;
            coeffVP += (1. - g)*sf.y;

            triplets_.push_back(SparseMatrix::Triplet(rowU, colU, coeffUV));
            triplets_.push_back(SparseMatrix::Triplet(rowV, colV, coeffUV));
            triplets_.push_back(SparseMatrix::Triplet(rowU, colP, coeffUP));
            triplets_.push_back(SparseMatrix::Triplet(rowV, colP, coeffVP));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Face &face = bd.face();
            const Vector2D& sf = bd.outwardNorm();
            const Vector2D& rf = bd.rFaceVec();

            const Scalar massFlux = rho_(face)*dot(u_(face), sf);
            const Scalar diffFlux = mu_(face)*dot(rf, sf)/dot(rf, rf);

            switch(u_.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                centralCoeffUV += diffFlux;

                rhs_(rowU) += (diffFlux - massFlux)*u_(bd.face()).x;
                rhs_(rowV) += (diffFlux - massFlux)*u_(bd.face()).y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeffUV += massFlux;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("CoupledEquation", "assembleMomentumEquation", "unrecognized or unspecified boundary type.");
            }

            switch(p_.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                rhs_(rowU) -= p_(bd.face())*sf.x;
                rhs_(rowV) -= p_(bd.face())*sf.y;
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT: case ScalarFiniteVolumeField::SYMMETRY:
                centralCoeffUP += sf.x;
                centralCoeffVP += sf.y;
                break;

            default:
                throw Exception("CoupledEquation", "assembleMomentumEquation", "unrecognized or unspecified boundary type.");
            }
        }

        d_(cell) = cell.volume()/centralCoeffUV;

        triplets_.push_back(SparseMatrix::Triplet(rowU, rowU, centralCoeffUV));
        triplets_.push_back(SparseMatrix::Triplet(rowV, rowV, centralCoeffUV));
        triplets_.push_back(SparseMatrix::Triplet(rowU, rowP, centralCoeffUP));
        triplets_.push_back(SparseMatrix::Triplet(rowV, rowP, centralCoeffVP));
    }

    interpolateFaces(fv::INVERSE_VOLUME, d_);
}

void CoupledEquation::assembleContinuityEquation()
{

    for(const Cell &cell: p_.grid.fluidCells())
    {
        const Index rowU = cell.globalIndex();
        const Index rowV = rowU + nActiveCells_;
        const Index rowP = rowV + nActiveCells_;
        Scalar centralCoeffP = 0.;
        Scalar centralCoeffU = 0.;
        Scalar centralCoeffV = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Index colU = nb.cell().globalIndex();
            const Index colV = colU + nActiveCells_;
            const Index colP = colV + nActiveCells_;
            Scalar coeffP = 0.;
            Scalar coeffU = 0.;
            Scalar coeffV = 0.;

            const Face &face = nb.face();
            const Vector2D& sf = nb.outwardNorm();
            const Vector2D& rc = nb.rCellVec();
            const Scalar df = d_.faces()[face.id()];

            //- Pressure poisson equation
            const Scalar diffFlux = df*dot(rc, sf)/dot(rc, rc);

            centralCoeffP -= diffFlux;
            coeffP += diffFlux;

            //- Source terms
            const Scalar g = nb.cell().volume()/(cell.volume() + nb.cell().volume());

            centralCoeffU -= g*sf.x;
            centralCoeffV -= g*sf.y;
            coeffU -= (1. - g)*sf.x;
            coeffV -= (1. - g)*sf.y;

            rhs_(rowP) += dot(g*d_(cell)*gradP_(cell) + (1. - g)*d_(nb.cell())*gradP_(nb.cell()), sf);

            triplets_.push_back(SparseMatrix::Triplet(rowP, colP, coeffP));
            triplets_.push_back(SparseMatrix::Triplet(rowP, colU, coeffU));
            triplets_.push_back(SparseMatrix::Triplet(rowP, colV, coeffV));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Face &face = bd.face();
            const Vector2D &sf = bd.outwardNorm();
            const Vector2D &rf = bd.rFaceVec();

            const Scalar massFlux = dot(u_(face), sf);
            const Scalar diffFlux = d_(face)*dot(rf, sf)/dot(rf, rf);

            switch(u_.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                rhs_(rowP) += massFlux;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeffU -= sf.x;
                centralCoeffV -= sf.y;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                break;
            default:
                throw Exception("CoupledEquation", "assembleContinuity", "unrecognized or unspecified boundary type.");
            }

            switch(p_.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                rhs_(rowP) -= diffFlux*p_(face);
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("CoupledEquation", "assembleContinuity", "unrecognized or unspecified boundary type.");
            }
        }

        triplets_.push_back(SparseMatrix::Triplet(rowP, rowP, centralCoeffP));
        triplets_.push_back(SparseMatrix::Triplet(rowP, rowU, centralCoeffU));
        triplets_.push_back(SparseMatrix::Triplet(rowP, rowV, centralCoeffV));
    }
}

void CoupledEquation::rhieChowInterpolation()
{
    computeGradient(fv::FACE_TO_CELL, p_, gradP_);

    for(const Face& face: u_.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        const Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());

        u_(face) = g*(u_(lCell) + d_(lCell)*gradP_(lCell)) + (1. - g)*(u_(rCell) + d_(rCell)*gradP_(rCell))
                - d_(face)*gradP_(face);
    }

    for(const Face& face: u_.grid.boundaryFaces())
    {
        switch(u_.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u_(face) = u_(face.lCell()) + d_(face.lCell())*gradP_(face.lCell())
                    - d_(face)*gradP_(face);
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            u_(face) = u_(face.lCell()) - dot(u_(face.lCell()), nWall)*nWall/nWall.magSqr();
        }
            break;

        default:
            throw Exception("CoupledEquation", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}
