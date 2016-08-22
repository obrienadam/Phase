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

    spMat_.setTolerance(1e-7);
    spMat_.setMaxIterations(1000);
    spMat_.setFill(0);
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
    phi_ = spMat_.solve(rhs_, phi_, SparseMatrix::IncompleteLUT);
    printf("Solved coupled Navier-Stokes equations. Error = %lf, number of iterations = %d.\n", spMat_.error(), spMat_.nIterations());

    for(const Cell &cell: u_.grid.activeCells())
        u_[cell.id()].x = phi_(cell.globalIndex());

    for(const Cell &cell: u_.grid.activeCells())
        u_[cell.id()].y = phi_(cell.globalIndex() + nActiveCells_);

    for(const Cell &cell: p_.grid.activeCells())
        p_[cell.id()] = phi_(cell.globalIndex() + 2*nActiveCells_);

    return spMat_.error();
}

//- Protected methods

void CoupledEquation::assembleMomentumEquation(Scalar timeStep)
{
    d_.fill(0.);

    for(const Cell &cell: u_.grid.fluidCells())
    {
        const size_t rowU = cell.globalIndex();
        const size_t rowV = rowU + nActiveCells_;
        const size_t rowP = rowV + nActiveCells_;
        Scalar centralCoeff = 0.;
        Scalar centralCoeffUP = 0.;
        Scalar centralCoeffVP = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const size_t colU = nb.cell().globalIndex();
            const size_t colV = colU + nActiveCells_;
            const size_t colP = colV + nActiveCells_;
            Scalar coeff = 0.;
            Scalar coeffUP = 0.;
            Scalar coeffVP = 0.;

            const Face &face = nb.face();
            const Vector2D& sf = nb.outwardNorm();
            const Vector2D& rc = nb.rCellVec();

            //- Time term
            centralCoeff += cell.volume()/timeStep;

            rhs_(rowU) += cell.volume()*u_.prev(0)[cell.id()].x/timeStep;
            rhs_(rowV) += cell.volume()*u_.prev(0)[cell.id()].y/timeStep;

            //- Advection term
            const Scalar massFlux = rho_.faces()[face.id()]*dot(u_.faces()[face.id()], sf);

            centralCoeff += std::max(massFlux, 0.);
            coeff += std::min(massFlux, 0.);

            //- Diffusion term
            const Scalar diffFlux = mu_.faces()[face.id()]*dot(rc, sf)/dot(rc, rc);

            centralCoeff += diffFlux;
            coeff -= diffFlux;

            //- Pressure term
            const Scalar g = nb.cell().volume()/(cell.volume() + nb.cell().volume());

            centralCoeffUP += g*sf.x;
            centralCoeffVP += g*sf.y;

            coeffUP += (1. - g)*sf.x;
            coeffVP += (1. - g)*sf.y;

            triplets_.push_back(SparseMatrix::Triplet(rowU, colU, coeff));
            triplets_.push_back(SparseMatrix::Triplet(rowV, colV, coeff));
            triplets_.push_back(SparseMatrix::Triplet(rowU, colP, coeffUP));
            triplets_.push_back(SparseMatrix::Triplet(rowV, colP, coeffVP));
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
                centralCoeff += diffFlux;

                rhs_(rowU) += (diffFlux - massFlux)*u_.faces()[bd.face().id()].x;
                rhs_(rowV) += (diffFlux - massFlux)*u_.faces()[bd.face().id()].y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += massFlux;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                break;

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
                centralCoeffUP += sf.x;
                centralCoeffVP += sf.y;
                break;

            default:
                throw Exception("CoupledEquation", "assembleMomentumEquation", "unrecognized or unspecified boundary type.");
            }
        }

        d_[cell.id()] = cell.volume()/centralCoeff;

        triplets_.push_back(SparseMatrix::Triplet(rowU, rowU, centralCoeff));
        triplets_.push_back(SparseMatrix::Triplet(rowV, rowV, centralCoeff));
        triplets_.push_back(SparseMatrix::Triplet(rowU, rowP, centralCoeffUP));
        triplets_.push_back(SparseMatrix::Triplet(rowV, rowP, centralCoeffVP));
    }

    interpolateFaces(fv::INVERSE_VOLUME, d_);
}

void CoupledEquation::assembleContinuityEquation()
{

    for(const Cell &cell: p_.grid.fluidCells())
    {
        const size_t rowU = cell.globalIndex();
        const size_t rowV = rowU + nActiveCells_;
        const size_t rowP = rowV + nActiveCells_;
        Scalar centralCoeff = 0.;
        Scalar centralCoeffU = 0.;
        Scalar centralCoeffV = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const size_t colU = nb.cell().globalIndex();
            const size_t colV = colU + nActiveCells_;
            const size_t colP = colV + nActiveCells_;
            Scalar coeff = 0.;
            Scalar coeffU = 0.;
            Scalar coeffV = 0.;

            const Face &face = nb.face();
            const Vector2D& sf = nb.outwardNorm();
            const Vector2D& rc = nb.rCellVec();
            const Scalar df = d_.faces()[face.id()];

            //- Pressure poisson equation
            const Scalar diffFlux = df*dot(rc, sf)/dot(rc, rc);

            centralCoeff -= diffFlux;
            coeff += diffFlux;

            //- Source terms
            const Scalar g = nb.cell().volume()/(cell.volume() + nb.cell().volume());

            centralCoeffU -= g*sf.x;
            centralCoeffV -= g*sf.y;
            coeffU -= (1. - g)*sf.x;
            coeffV -= (1. - g)*sf.y;

            rhs_(rowP) += dot(g*d_[cell.id()]*gradP_[cell.id()] + (1. - g)*d_[nb.cell().id()]*gradP_[nb.cell().id()], sf);

            triplets_.push_back(SparseMatrix::Triplet(rowP, colP, coeff));
            triplets_.push_back(SparseMatrix::Triplet(rowP, colU, coeffU));
            triplets_.push_back(SparseMatrix::Triplet(rowP, colV, coeffV));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D &sf = bd.outwardNorm();
            const Scalar massFlux = dot(u_.faces()[bd.face().id()], sf);
            const Scalar diffFlux = d_.faces()[bd.face().id()]*dot(bd.rFaceVec(), sf)/dot(bd.rFaceVec(), bd.rFaceVec());

            switch(u_.boundaryType(bd.face().id()))
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

            switch(p_.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                rhs_(rowP) -= diffFlux*p_.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:

                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("CoupledEquation", "assembleContinuity", "unrecognized or unspecified boundary type.");
            }
        }

        triplets_.push_back(SparseMatrix::Triplet(rowP, rowP, centralCoeff));
        triplets_.push_back(SparseMatrix::Triplet(rowP, rowU, centralCoeffU));
        triplets_.push_back(SparseMatrix::Triplet(rowP, rowV, centralCoeffV));
    }
}

void CoupledEquation::rhieChowInterpolation()
{
    computeGradient(fv::GREEN_GAUSS_CELL_CENTERED, p_, gradP_);

    for(const Face& face: u_.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        const Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());

        const Vector2D sf = face.outwardNorm(lCell.centroid());
        const Vector2D rc = rCell.centroid() - lCell.centroid();

        u_.faces()[face.id()] -= d_.faces()[face.id()]*((p_[rCell.id()] - p_[lCell.id()])*sf/dot(rc, rc) - gradP_.faces()[face.id()]);

        u_.faces()[face.id()] = g*u_[lCell.id()] + (1. - g)*u_[rCell.id()]
                - d_.faces()[face.id()]*(p_[rCell.id()] - p_[lCell.id()])*rc/dot(rc, rc)
                + g*d_[lCell.id()]*gradP_[lCell.id()] + (1. - g)*d_[rCell.id()]*gradP_[rCell.id()];
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
