#include "SurfaceTensionForce.h"
#include "GhostCellStencil.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input,
                                         const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                         const std::weak_ptr<ImmersedBoundary> &ib)
        :
        VectorFiniteVolumeField(grid, "ft", Vector2D(0., 0.), true, false),
        ib_(ib),
        kappa_(std::make_shared<ScalarFiniteVolumeField>(grid_, "kappa")),
        gammaTilde_(std::make_shared<ScalarFiniteVolumeField>(grid_, "gammaTilde")),
        gradGammaTilde_(std::make_shared<ScalarGradient>(*gammaTilde_)),
        n_(std::make_shared<VectorFiniteVolumeField>(grid_, "n"))
{
    //- Input properties
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");

    //- Determine which patches contact angles will be enforced on
    for (const auto &input: input.boundaryInput().get_child("Boundaries.gamma"))
    {
        if (input.first == "*" || !grid_->hasPatch(input.first))
            continue;

        patchContactAngles_.insert(std::make_pair(
                grid_->patch(input.first).id(),
                input.second.get<Scalar>("contactAngle", 90) * M_PI / 180.
        ));
    }

    //- Determine which IBs contact angles will be enforced on
    if (ib_.lock())
        for (const auto &ibObj: *ib_.lock())
        {
            ibContactAngles_[ibObj->id()] = input.boundaryInput().get<Scalar>(
                    "ImmersedBoundaries." + ibObj->name() + ".gamma.contactAngle", 90) * M_PI / 180.;
        }
}

void SurfaceTensionForce::computeInterfaceNormals()
{
    const VectorFiniteVolumeField &gradGammaTilde = *gradGammaTilde_;
    VectorFiniteVolumeField &n = *n_;

    for (const Cell &cell: grid_->cellZone("fluid"))
        n(cell) = gradGammaTilde(cell).magSqr() >= eps_ * eps_ ? -gradGammaTilde(cell).unitVec() : Vector2D(0., 0.);

    //- Boundary faces set from contact line orientation
    for (const Patch &patch: grid_->patches())
    {
        for (const Face &face: patch)
        {
            if (n(face.lCell()).magSqr() == 0.)
            {
                n(face) = Vector2D(0., 0.);
                continue;
            }

            Vector2D ns = -face.outwardNorm(face.lCell().centroid());
            Vector2D ts = (n(face.lCell()) - dot(n(face.lCell()), ns) * ns).unitVec();

            Scalar theta = this->theta(patch);

            n(face) = ns * std::cos(theta) + ts * std::sin(theta);
        }
    }

    grid_->sendMessages(n);
}

Scalar SurfaceTensionForce::theta(const Patch &patch) const
{
    auto it = patchContactAngles_.find(patch.id());
    return it != patchContactAngles_.end() ? it->second : M_PI_2;
}

Scalar SurfaceTensionForce::theta(const ImmersedBoundaryObject &ibObj) const
{
    auto it = ibContactAngles_.find(ibObj.id());
    return it != ibContactAngles_.end() ? it->second : M_PI_2;
}

Vector2D
SurfaceTensionForce::contactLineNormal(const Cell &cell, const Point2D &pt, const ImmersedBoundaryObject &ibObj) const
{
    const Vector2D &n = (*n_)(cell);
    Vector2D ns = -ibObj.nearestEdgeNormal(pt);
    Vector2D ts = (n - dot(n, ns) * ns).unitVec();

    if (n.magSqr() == 0)
        return n;

    Scalar theta = this->theta(ibObj);
    return ns * std::cos(theta) + ts * std::sin(theta);
}

void SurfaceTensionForce::smoothGammaField(const ScalarFiniteVolumeField &gamma)
{
    gammaTilde_->fill(0);

    if (ib_.lock())
    {
        CellGroup cellsToSmooth = grid_->localActiveCells() - ib_.lock()->solidCells();

        smooth(gamma,
               cellsToSmooth,
               grid_->globalCellGroup(cellsToSmooth),
               kernelWidth_,
               *gammaTilde_,
               [](const Cell &cell, const Cell &kCell, Scalar e)
               {
                   Scalar r = (cell.centroid() - kCell.centroid()).mag() / e;
                   return r < 1. ? std::cos(M_PI * r) + 1. : 0.;
               });

        grid_->sendMessages(*gammaTilde_);
    }
    else
    {
        smooth(gamma,
               grid_->localActiveCells(),
               grid_->globalActiveCells(),
               kernelWidth_,
               *gammaTilde_,
               [](const Cell &cell, const Cell &kCell, Scalar e)
               {
                   Scalar r = (cell.centroid() - kCell.centroid()).mag() / e;
                   return r < 1. ? std::cos(M_PI * r) + 1. : 0.;
               });
    }

    gammaTilde_->setBoundaryFaces();
}

Equation<Scalar> SurfaceTensionForce::contactLineBcs(ScalarFiniteVolumeField &gamma)
{
    Equation<Scalar> eqn(gamma);

    for (auto ibObj: *ib_.lock())
    {
        switch (ibObj->type())
        {
            case ImmersedBoundaryObject::GHOST_CELL:
                eqn += ibObj->contactLineBcs(gamma, theta(*ibObj));
                break;
            case ImmersedBoundaryObject::QUADRATIC:
            {
                Scalar theta = this->theta(*ibObj);

                for (const Cell &cell: ibObj->ibCells())
                {
                    Vector2D wn = -ibObj->nearestEdgeNormal(cell.centroid());

                    Ray2D r1 = Ray2D(cell.centroid(), wn.rotate(M_PI_2 - theta));
                    Ray2D r2 = Ray2D(cell.centroid(), wn.rotate(theta - M_PI_2));

                    GhostCellStencil m1(cell, ibObj->shape().intersections(r1)[0], r1.r(), *grid_);
                    GhostCellStencil m2(cell, ibObj->shape().intersections(r2)[0], r2.r(), *grid_);
                    Scalar g1 = m1.bpValue(gamma);
                    Scalar g2 = m2.bpValue(gamma);

                    if (std::abs(g1 - g2) > 1e-8)
                    {
                        if (g2 < g1)
                            std::swap(m1, m2);
                    }
                    else
                    {
                        Vector2D grad1 = m1.bpGrad(gamma);
                        Vector2D grad2 = m2.bpGrad(gamma);

                        if (g2 + dot(grad2, r2.r()) < g1 + dot(grad1, r1.r()))
                            std::swap(m1, m2);
                    }

                    if (theta > M_PI_2)
                        eqn.add(m1.cell(), m1.neumannCells(), m1.neumannCoeffs());
                    else
                        eqn.add(m2.cell(), m2.neumannCells(), m2.neumannCoeffs());
                }

                for (const Cell &cell: ibObj->solidCells())
                    eqn.add(cell, cell, 1.);
            }

                break;

            case ImmersedBoundaryObject::HIGH_ORDER:
                eqn += ibObj->contactLineBcs(gamma, theta(*ibObj));

                break;
            default:
                throw Exception("SurfaceTensionForce", "contactLineBcs", "unrecognized immersed boundary object type.");
        }
    }

    return eqn;
}
