#include <ImmersedBoundary/QuadraticImmersedBoundaryObject.h>
#include <ImmersedBoundary/HighOrderImmersedBoundaryObject.h>
#include "Celeste.h"
#include "Algorithm.h"
#include "GhostCellImmersedBoundaryObject.h"

Celeste::Celeste(const Input &input,
                 const ImmersedBoundary &ib,
                 ScalarFiniteVolumeField &gamma,
                 const ScalarFiniteVolumeField &rho,
                 const ScalarFiniteVolumeField &mu,
                 const VectorFiniteVolumeField &u,
                 const ScalarGradient &gradGamma)
        :
        SurfaceTensionForce(input, ib, gamma, rho, mu, u, gradGamma)
{
    constructMatrices();
}

void Celeste::computeFaces()
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    auto &ft = *this;
    auto &kappa = *kappa_;

    for (const Face &face: gamma_.grid()->faces())
        ft(face) = sigma_ * kappa(face) * gradGamma_(face);
}

void Celeste::computeFaces(const ImmersedBoundary &ib)
{
    computeGradGammaTilde(ib);
    computeInterfaceNormals();
    computeCurvature(ib);

    auto &ft = *this;
    auto &kappa = *kappa_;

    for (const Face &face: gamma_.grid()->faces())
        ft(face) = sigma_ * kappa(face) * gradGamma_(face);
}

void Celeste::compute()
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    auto &ft = *this;
    auto &kappa = *kappa_;

    ft.fill(Vector2D(0., 0.));
    for (const Cell &cell: gamma_.grid()->cellZone("fluid"))
        ft(cell) = sigma_ * kappa(cell) * gradGamma_(cell);

    ft.interpolateFaces();
}

void Celeste::compute(const ImmersedBoundary &ib)
{
    computeFaces(ib);

    auto &ft = *this;
    const auto &kappa = *kappa_;

    for (const Cell &cell: grid_->cellZone("fluid"))
        ft(cell) = sigma_ * kappa(cell) * gradGamma_(cell);
}

void Celeste::constructMatrices()
{
    kappaStencils_.resize(grid_->cells().size());
    gradGammaTildeStencils_.resize(grid_->cells().size());

    for (const Cell &cell: grid_->localActiveCells())
    {
        kappaStencils_[cell.id()] = CelesteStencil(cell, false);
        gradGammaTildeStencils_[cell.id()] = CelesteStencil(cell, true);
    }
}

Equation<Scalar> Celeste::contactLineBcs(const ImmersedBoundary &ib)
{
    Equation<Scalar> eqn(gamma_);

    for (auto ibObj: ib)
    {
        switch (ibObj->type())
        {
            case ImmersedBoundaryObject::GHOST_CELL:
                eqn += ibObj->contactLineBcs(gamma_, getTheta(*ibObj));
                break;
            case ImmersedBoundaryObject::QUADRATIC:
                for (const Cell &cell: ibObj->ibCells())
                {
                    Vector2D wn = -ibObj->nearestEdgeNormal(ibObj->nearestIntersect(cell.centroid()));

                    Scalar theta = getTheta(*ibObj);

                    Ray2D r1 = Ray2D(cell.centroid(), wn.rotate(M_PI_2 - theta));
                    Ray2D r2 = Ray2D(cell.centroid(), wn.rotate(theta - M_PI_2));

                    GhostCellStencil m1(cell, ibObj->shape().intersections(r1)[0], r1.r(), *grid());
                    GhostCellStencil m2(cell, ibObj->shape().intersections(r2)[0], r2.r(), *grid());

                    if (theta < M_PI_2)
                    {
                        if (m1.ipValue(gamma_) > m2.ipValue(gamma_))
                            eqn.add(m1.cell(), m1.neumannCells(), m1.neumannCoeffs());
                        else
                            eqn.add(m2.cell(), m2.neumannCells(), m2.neumannCoeffs());
                    }
                    else
                    {
                        if (m1.ipValue(gamma_) < m2.ipValue(gamma_))
                            eqn.add(m1.cell(), m1.neumannCells(), m1.neumannCoeffs());
                        else
                            eqn.add(m2.cell(), m2.neumannCells(), m2.neumannCoeffs());
                    }
                }

                for (const Cell &cell: ibObj->solidCells())
                    eqn.add(cell, cell, 1.);

                break;

            case ImmersedBoundaryObject::HIGH_ORDER:
                for(auto ibObj: ib)
                    eqn += ibObj->contactLineBcs(gamma_, getTheta(*ibObj));

                break;
            default:
                throw Exception("Celeste", "contactLineBcs", "unrecognized immersed boundary object type.");
        }
    }

    return eqn;
}

//- Protected methods

void Celeste::computeGradGammaTilde()
{
    smoothGammaField();

    auto &gammaTilde = *gammaTilde_;
    auto &gradGammaTilde = *gradGammaTilde_;

    gradGammaTilde_->fill(Vector2D(0., 0.));
    for (const Cell &cell: gradGamma_.grid()->cellZone("fluid"))
        gradGammaTilde(cell) = gradGammaTildeStencils_[cell.id()].grad(gammaTilde);
}

void Celeste::computeGradGammaTilde(const ImmersedBoundary &ib)
{
    smoothGammaField(ib);

    auto &gammaTilde = *gammaTilde_;
    auto &gradGammaTilde = *gradGammaTilde_;

    gradGammaTilde_->fill(Vector2D(0., 0.));
    for (const Cell &cell: gradGamma_.grid()->cellZone("fluid"))
        gradGammaTilde(cell) = gradGammaTildeStencils_[cell.id()].grad(gammaTilde);
}

void Celeste::computeCurvature()
{
    auto &n = *n_;
    auto &kappa = *kappa_;

    for (const Cell &cell: grid_->cellZone("fluid"))
        kappa(cell) = kappaStencils_[cell.id()].div(n);

    grid_->sendMessages(kappa);
    kappa.interpolateFaces();
}

void Celeste::computeCurvature(const ImmersedBoundary &ib)
{
    updateStencils(ib);
    kappa_->fill(0.);

    auto &n = *n_;
    auto &kappa = *kappa_;
    const auto &gradGammaTilde = *gradGammaTilde_;

    for (const Cell &cell: grid_->cellZone("fluid"))
        if (gradGammaTilde(cell).magSqr() > 0.)
            kappa(cell) = kappaStencils_[cell.id()].kappa(n, ib, *this);

    grid_->sendMessages(kappa);
    kappa.interpolateFaces();

    for(const Face& face: grid_->interiorFaces())
    {
        if(ib.ibObj(face.lCell().centroid()))
        {
            kappa(face) = kappa(face.rCell());
        }
        else if(ib.ibObj(face.rCell().centroid()))
        {
            kappa(face) = kappa(face.lCell());
        }
    }
}

void Celeste::updateStencils(const ImmersedBoundary &ib)
{
//    auto updateRequired = [&ib](const CelesteStencil &st) {
//        if (st.truncated())
//            return true;
//
//        for (const InteriorLink &nb: st.cell().neighbours())
//            if (ib.ibObj(nb.cell().centroid()))
//                return true;
//
//        for (const CellLink &dg: st.cell().diagonals())
//            if (ib.ibObj(dg.cell().centroid()))
//                return true;
//
//        return false;
//    };

    for (const Cell &cell: grid_->cellZone("fluid"))
    {
        CelesteStencil &st = kappaStencils_[cell.id()];

    //    if (updateRequired(st))
            st.init(ib);
    }
}
