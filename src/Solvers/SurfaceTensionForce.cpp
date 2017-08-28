#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input,
                                         const ImmersedBoundary &ib,
                                         const ScalarFiniteVolumeField &gamma,
                                         const ScalarFiniteVolumeField &rho,
                                         const ScalarFiniteVolumeField &mu,
                                         const VectorFiniteVolumeField &u,
                                         const ScalarGradient &gradGamma)
        :
        VectorFiniteVolumeField(gamma.gridPtr(), "ft", Vector2D(0., 0.), true, false),
        ib_(ib),
        gamma_(gamma),
        rho_(rho),
        mu_(mu),
        u_(u),
        gradGamma_(gradGamma),
        kappa_(std::make_shared<ScalarFiniteVolumeField>(grid_, "kappa")),
        gammaTilde_(std::make_shared<ScalarFiniteVolumeField>(grid_, "gammaTilde")),
        gradGammaTilde_(std::make_shared<ScalarGradient>(*gammaTilde_)),
        n_(std::make_shared<VectorFiniteVolumeField>(grid_, "n"))
{
    //- Input properties
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");

    //- Determine which patches contact angles will be enforced on
    for (const auto &input: input.boundaryInput().get_child("Boundaries." + gamma.name()))
    {
        if (input.first == "*" || !grid_->hasPatch(input.first))
            continue;

        patchContactAngles_.insert(std::make_pair(
                grid_->patch(input.first).id(),
                input.second.get<Scalar>("contactAngle", 90) * M_PI / 180.
        ));
    }

    //- Determine which IBs contact angles will be enforced on

    if (input.boundaryInput().find("ImmersedBoundaries") == input.boundaryInput().not_found())
        return;

    for (const auto &input: input.boundaryInput().get_child("ImmersedBoundaries"))
    {
        ibContactAngles_.insert(std::make_pair(
                ib_.ibObj(input.first).id(),
                input.second.get<Scalar>(gamma.name() + ".contactAngle", 90) * M_PI / 180.
        ));
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
            Scalar theta = getTheta(patch);

            if (n(face.lCell()) == Vector2D(0., 0.))
                n(face) = Vector2D(0., 0.);
            else
            {
                Vector2D t = face.norm().tangentVec().unitVec();
                t = dot(t, n(face.lCell())) > 0. ? t : -t;
                Scalar phi = cross(t, face.outwardNorm()) < 0. ? M_PI_2 - theta : theta - M_PI_2;
                n(face) = t.rotate(phi);
            }
        }
    }

    grid_->sendMessages(n);
}

Vector2D SurfaceTensionForce::contactLineNormal(const Cell &lCell,
                                                const Cell &rCell,
                                                const ImmersedBoundaryObject &ibObj)
{
    LineSegment2D ln = ibObj.intersectionLine(LineSegment2D(lCell.centroid(), rCell.centroid()));
    Vector2D en = ibObj.nearestEdgeNormal(ln.ptB());

    const VectorFiniteVolumeField &n = *n_;

    if(n(lCell) == Vector2D(0., 0.))
        return Vector2D(0., 0.);

    Vector2D t = en.tangentVec().unitVec();
    Scalar theta = getTheta(ibObj);

    t = dot(t, n(lCell)) > 0. ? t : -t;
    Scalar phi = cross(t, en) < 0. ? M_PI_2 - theta : theta - M_PI_2;
    return t.rotate(phi);
}

void SurfaceTensionForce::smoothGammaField()
{
//    smooth(gamma_, grid().localActiveCells(), kernelWidth_, *gammaTilde_,
//           [](const Cell& cell, const Cell& kCell, Scalar e) {
//               Vector2D r = (cell.centroid() - kCell.centroid()).abs() / e;
//               Scalar dx = r.x < 1. ? pow(r.x, 4) - 2*pow(r.x, 2) + 1. : 0.;
//               Scalar dy = r.y < 1. ? pow(r.y, 4) - 2*pow(r.y, 2) + 1. : 0.;
//               return dx * dy;
//           });

//    smooth(gamma_,
//           grid().localActiveCells(),
//           grid().globalActiveCells(),
//           kernelWidth_,
//           *gammaTilde_,
//           [](const Cell& cell, const Cell& kCell, Scalar e) {
//               Scalar r = (cell.centroid() - kCell.centroid()).mag() / e;
//               return r < 1. ? pow(r, 4) - 2 * pow(r, 2) + 1 : 0.;
//           });

//    smooth(gamma_, grid().localActiveCells(), kernelWidth_, *gammaTilde_,
//    [](const Cell& cell, const Cell &kCell, Scalar e) {
//        Scalar r = (cell.centroid() - kCell.centroid()).mag() / e;
//        return r < 1. ? pow(1. - r*r, 3) : 0.;
//    });
    //- This kernel is better
    smooth(gamma_, grid().localActiveCells(), grid().globalActiveCells(), kernelWidth_, *gammaTilde_,
           [](const Cell &cell, const Cell &kCell, Scalar e) {
               Scalar r = (cell.centroid() - kCell.centroid()).mag() / e;
               return r < 1. ? std::cos(M_PI * r) + 1. : 0.;
           });

    grid_->sendMessages(*gammaTilde_);
    gammaTilde_->setBoundaryFaces();
}