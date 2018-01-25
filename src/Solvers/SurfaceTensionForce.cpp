#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input,
                                         const ImmersedBoundary &ib,
                                         ScalarFiniteVolumeField &gamma,
                                         const ScalarFiniteVolumeField &rho,
                                         const ScalarFiniteVolumeField &mu,
                                         const VectorFiniteVolumeField &u,
                                         const ScalarGradient &gradGamma)
        :
        VectorFiniteVolumeField(gamma.grid(), "ft", Vector2D(0., 0.), true, false),
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
            Vector2D t = face.norm().tangentVec().unitVec();
            t = dot(t, n(face.lCell())) > 0. ? t : -t;
            Scalar phi = cross(t, face.outwardNorm()) < 0. ? M_PI_2 - theta : theta - M_PI_2;
            n(face) = t.rotate(phi) * n(face.lCell()).magSqr();
        }
    }

    grid_->sendMessages(n);
}

Vector2D SurfaceTensionForce::contactLineNormal(const Cell &lCell,
                                                const Cell &rCell,
                                                const ImmersedBoundaryObject &ibObj) const
{
    //LineSegment2D ln = ibObj.intersectionLine(LineSegment2D(lCell.centroid(), rCell.centroid()));
    Vector2D en = ibObj.nearestEdgeNormal(rCell.centroid());

    const VectorFiniteVolumeField &n = *n_;

    if (n(lCell) == Vector2D(0., 0.))
        return Vector2D(0., 0.);

    Vector2D t = en.tangentVec().unitVec();
    Scalar theta = getTheta(ibObj);

    t = dot(t, n(lCell)) > 0. ? t : -t;
    Scalar phi = cross(t, en) < 0. ? M_PI_2 - theta : theta - M_PI_2;
    return t.rotate(phi) * n(lCell).magSqr();
}

Vector2D
SurfaceTensionForce::contactLineNormal(const Cell &cell, const Point2D &pt, const ImmersedBoundaryObject &ibObj) const
{
    const Vector2D &n = (*n_)(cell);
    Vector2D ns = -ibObj.nearestEdgeNormal(pt);
    Vector2D ts = (n - dot(n, ns) * ns).unitVec();

    if(n.magSqr() == 0)
        return n;

    Scalar theta = getTheta(ibObj);
    return ns * std::cos(theta) + ts * std::sin(theta);
}

Vector2D SurfaceTensionForce::contactLineNormal(const Cell &cell, const ImmersedBoundaryObject &ibObj) const
{
    //LineSegment2D ln = ibObj.intersectionLine(LineSegment2D(lCell.centroid(), rCell.centroid()));
    Vector2D en = ibObj.nearestEdgeNormal(cell.centroid());

    const VectorFiniteVolumeField &n = *n_;

    Vector2D t = en.tangentVec().unitVec();
    Scalar theta = getTheta(ibObj);

    t = dot(t, n(cell)) > 0. ? t : -t;
    Scalar phi = cross(t, en) < 0. ? M_PI_2 - theta : theta - M_PI_2;
    return t.rotate(phi) * n(cell).magSqr();
}

Vector2D SurfaceTensionForce::contactLineNormal(const Cell &cell, const ImmersedBoundary &ib) const
{
    auto ibObj = ib.nearestIbObj(cell.centroid());
    return contactLineNormal(cell, *ibObj);
}

void SurfaceTensionForce::smoothGammaField()
{
    const CellGroup &cellsToSmooth = grid_->globalActiveCells();

    smooth(gamma_, cellsToSmooth,
           cellsToSmooth,
           kernelWidth_,
           *gammaTilde_,
           [](const Cell &cell, const Cell &kCell, Scalar e) {
               Scalar r = (cell.centroid() - kCell.centroid()).mag() / e;
               return r < 1. ? std::cos(M_PI * r) + 1. : 0.;
           });

    grid_->sendMessages(*gammaTilde_);
    gammaTilde_->setBoundaryFaces();
}

void SurfaceTensionForce::smoothGammaField(const ImmersedBoundary &ib)
{
    CellGroup cellsToSmooth = grid_->localActiveCells() - ib.solidCells();

    smooth(gamma_, cellsToSmooth,
           grid_->globalCellGroup(cellsToSmooth),
           kernelWidth_,
           *gammaTilde_,
           [](const Cell &cell, const Cell &kCell, Scalar e) {
               Scalar r = (cell.centroid() - kCell.centroid()).mag() / e;
               return r < 1. ? std::cos(M_PI * r) + 1. : 0.;
           });

    grid_->sendMessages(*gammaTilde_);
    gammaTilde_->setBoundaryFaces();
}
