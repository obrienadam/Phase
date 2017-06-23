#include "ContinuumSurfaceForce.h"
#include "Interpolation.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

ContinuumSurfaceForce::ContinuumSurfaceForce(const Input &input,
                                             const ScalarFiniteVolumeField& gamma,
                                             const ScalarFiniteVolumeField& rho,
                                             const ScalarFiniteVolumeField& mu,
                                             const VectorFiniteVolumeField& u,
                                             VectorFiniteVolumeField& gradGamma)
        :
        SurfaceTensionForce(input, gamma, rho, mu, u, gradGamma),
        gammaTilde_(std::make_shared<ScalarFiniteVolumeField>(grid_, "gammaTilde")),
        gradGammaTilde_(std::make_shared<VectorFiniteVolumeField>(grid_, "gradGammaTilde"))
{
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");
    curvatureCutoffTolerance_ = input.caseInput().get<Scalar>("Solver.curvatureCutoffTolerance", 1e-10);

    gammaTilde_->copyBoundaryTypes(gamma);
    gammaTilde_->setBoundaryFaces(ScalarFiniteVolumeField::FIXED, [this](const Face& face){
        return gamma_(face);
    });

    constructSmoothingKernels();
}

void ContinuumSurfaceForce::compute()
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    VectorFiniteVolumeField &ft = *this;
    const ScalarFiniteVolumeField &kappa = *kappa_;

    ft.fill(Vector2D(0., 0.));

    for(const Cell &cell: grid_->cellZone("fluid"))
        ft(cell) = sigma_*kappa(cell)*gradGamma_(cell);

    for(const Face &face: gamma_.grid().faces())
        ft(face) = sigma_*kappa(face)*gradGamma_(face);

    grid_->sendMessages(ft);
}

void ContinuumSurfaceForce::constructSmoothingKernels()
{
    cellRangeSearch_ = gamma_.grid().constructSmoothingKernels(kernelWidth_);
}

void ContinuumSurfaceForce::registerSubFields(Solver& solver)
{
    solver.registerField(gammaTilde_);
    solver.registerField(gradGammaTilde_);
    solver.registerField(n_);
    solver.registerField(kappa_);
}

//- Private methods

void ContinuumSurfaceForce::computeGradGammaTilde()
{
    *gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    grid_->sendMessages(*gammaTilde_);
    computeGradient(fv::FACE_TO_CELL, grid_->cellZone("fluid"), *gammaTilde_, *gradGammaTilde_);
}

void ContinuumSurfaceForce::computeInterfaceNormals()
{
    VectorFiniteVolumeField &n = *n_;
    const VectorFiniteVolumeField &gradGammaTilde = *gradGammaTilde_;

    for (const Cell &cell: grid_->cells())
    {
        Vector2D m = -gradGammaTilde(cell);
        n(cell) = m.magSqr() > curvatureCutoffTolerance_ ? m.unitVec(): Vector2D(0., 0.);
    }

    grid_->sendMessages(n);

    for(const Patch& patch: grid_->patches())
    {
        if(isContactLinePatch(patch))
        {
            for(const Face& face: patch)
            {
                Vector2D sf = face.outwardNorm(face.lCell().centroid());
                n(face) = computeContactLineNormal(gradGammaTilde(face.lCell()), sf, u_(face.lCell()), theta());
            }
        }
        else
        {
            for(const Face& face: patch)
                n(face) = n(face.lCell());
        }
    }
}

void ContinuumSurfaceForce::computeCurvature()
{
    ScalarFiniteVolumeField &kappa = *kappa_;
    const VectorFiniteVolumeField &n = *n_;

    for(const Cell &cell: grid_->cellZone("fluid"))
    {
        Scalar &k = kappa(cell) = 0.;

        for(const InteriorLink &nb: cell.neighbours())
            k += dot(n(nb.face()), nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
            k += dot(n(bd.face()), bd.outwardNorm());

        k /= cell.volume();
    }

    grid_->sendMessages(kappa);
    interpolateCurvatureFaces();
}

void ContinuumSurfaceForce::interpolateCurvatureFaces()
{
    ScalarFiniteVolumeField &kappa = *kappa_;

    for(const Face &face: grid_->interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();
        const Vector2D rc = rCell.centroid() - lCell.centroid();

        const Scalar gradGammaMagSqr = ((gamma_(rCell) - gamma_(lCell))*rc/dot(rc, rc)).magSqr();

        if(gradGammaMagSqr < curvatureCutoffTolerance_)
            kappa(face) = 0.;
        else
        {
            Scalar g = rCell.volume()/(rCell.volume() + lCell.volume());
            kappa(face) = g*kappa(lCell) + (1. - g)*kappa(rCell);
        }
    }

    kappa.setBoundaryFaces();
}
