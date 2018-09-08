#include "FiniteVolume/ImmersedBoundary/ImmersedBoundary.h"
#include "FiniteVolumeGrid2D/BilinearInterpolator.h"

#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input,
                                         const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                         const std::shared_ptr<const CellGroup> &fluid)
    :
      grid_(grid),
      fluid_(fluid),
      fst_(std::make_shared<VectorFiniteVolumeField>(grid, "fst", 0., true, false, fluid)),
      kappa_(std::make_shared<ScalarFiniteVolumeField>(grid, "kappa", 0., true, false, fluid)),
      gammaTilde_(std::make_shared<ScalarFiniteVolumeField>(grid, "gammaTilde", 0., true, false, fluid)),
      gradGammaTilde_(std::make_shared<ScalarGradient>(*gammaTilde_, fluid)),
      n_(std::make_shared<VectorFiniteVolumeField>(grid, "n", Vector2D(0., 0.), true, false, fluid))
{
    //- Input properties
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");

    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");

    for(const Cell &cell: *fluid_)
        kernels_.push_back(SmoothingKernel(cell, kernelWidth_));

    //- Determine which patches contact angles will be enforced on
    for (const FaceGroup &patch: grid->patches())
        patchContactAngles_[patch.name()] = input.boundaryInput().get<Scalar>(
                    "Boundaries.gamma." + patch.name() + ".contactAngle", 90.) * M_PI / 180.;
}

void SurfaceTensionForce::setCellGroup(const std::shared_ptr<const CellGroup> &fluid)
{
    fluid_ = fluid;
    kappa_->setCellGroup(fluid);
    gammaTilde_->setCellGroup(fluid);
    n_->setCellGroup(fluid);
    gradGammaTilde_->setCellGroup(fluid);
}

void SurfaceTensionForce::computeInterfaceNormals()
{
    const VectorFiniteVolumeField &gradGammaTilde = *gradGammaTilde_;
    VectorFiniteVolumeField &n = *n_;

    for (const Cell &cell: *fluid_)
        n(cell) = gradGammaTilde(cell).magSqr() >= eps_ * eps_ ? -gradGammaTilde(cell).unitVec() : Vector2D(0., 0.);

    //- Boundary faces set from contact line orientation
    for (const FaceGroup &patch: n.grid()->patches())
    {
        Scalar theta = this->theta(patch);

        for (const Face &face: patch)
        {
            if (n(face.lCell()).magSqr() == 0.)
            {
                n(face) = Vector2D(0., 0.);
                continue;
            }

            Vector2D ns = -face.outwardNorm(face.lCell().centroid()).unitVec();
            Vector2D ts = (n(face.lCell()) - dot(n(face.lCell()), ns) * ns).unitVec();

            n(face) = ns * std::cos(theta) + ts * std::sin(theta);
        }
    }

    grid_->sendMessages(n);
}

Scalar SurfaceTensionForce::theta(const FaceGroup &patch) const
{
    auto it = patchContactAngles_.find(patch.name());
    return it != patchContactAngles_.end() ? it->second : M_PI_2;
}

void SurfaceTensionForce::smoothGammaField(const ScalarFiniteVolumeField &gamma)
{
    auto &gammaTilde = *gammaTilde_;

    gammaTilde.fill(0.);
    for(const auto &k: kernels_)
        gammaTilde(k.cell()) = k.eval(gamma);

    grid_->sendMessages(gammaTilde);
    gammaTilde.setBoundaryFaces();
}

//Vector2D SurfaceTensionForce::computeCapillaryForce(const ScalarFiniteVolumeField &gamma,
//                                                    const ImmersedBoundaryObject &ibObj) const
//{
//    typedef std::tuple<Point2D, Scalar> IbPoint;

//    std::vector<IbPoint> ibPoints(ibObj.ibCells().size());
//    std::transform(ibObj.ibCells().begin(), ibObj.ibCells().end(), ibPoints.begin(),
//                   [this, &gamma, &ibObj](const Cell &cell)
//    {
//        Point2D pt = ibObj.nearestIntersect(cell.centroid());
//        Scalar g = BilinearInterpolator(grid_, pt)(gamma);
//        return std::make_tuple(pt, g);
//    });

//    ibPoints = grid_->comm().gatherv(grid_->comm().mainProcNo(), ibPoints);

//    Vector2D force = Vector2D(0., 0.);

//    if (grid_->comm().isMainProc())
//    {
//        std::sort(ibPoints.begin(), ibPoints.end(), [&ibObj](const IbPoint &lhs, const IbPoint &rhs)
//        {
//            return (std::get<0>(lhs) - ibObj.shape().centroid()).angle()
//                    < (std::get<0>(rhs) - ibObj.shape().centroid()).angle();
//        });

//        for (int i = 0; i < ibPoints.size(); ++i)
//        {
//            const auto &a = ibPoints[i];
//            const auto &b = ibPoints[(i + 1) % ibPoints.size()];

//            if ((std::get<1>(a) < 0.5) != (std::get<1>(b) <= 0.5))
//            {
//                Scalar alpha = (0.5 - std::get<1>(b)) / (std::get<1>(a) - std::get<1>(b));
//                Point2D xcl = ibObj.nearestIntersect(alpha * std::get<0>(a) + (1. - alpha) * std::get<0>(b));
//                Vector2D tcl = contactLineTangent(
//                            std::get<1>(a) < std::get<1>(b) ?
//                                std::get<0>(a) - std::get<0>(b) : std::get<0>(b) - std::get<0>(a),
//                            xcl,
//                            ibObj);

//                force += sigma_ * tcl;
//            }
//        }

//        std::cout << "Capillary force = " << force << "\n";
//    }

//    return grid_->comm().broadcast(grid_->comm().mainProcNo(), force);
//}

//FiniteVolumeEquation<Scalar> SurfaceTensionForce::contactLineBcs(ScalarFiniteVolumeField &gamma)
//{
//    FiniteVolumeEquation<Scalar> eqn(gamma);

//    for (auto ibObj: *ib_.lock())
//    {
//        switch (ibObj->type())
//        {
//            case ImmersedBoundaryObject::GHOST_CELL:
//            case ImmersedBoundaryObject::QUADRATIC:
//            {
//                Scalar theta = this->theta(*ibObj);

//                for (const Cell &cell: ibObj->ibCells())
//                {
//                    Vector2D nw = -ibObj->nearestEdgeNormal(cell.centroid()).unitVec();
//                    Vector2D cl1 = nw.rotate(M_PI_2 - theta);
//                    Vector2D cl2 = nw.rotate(theta - M_PI_2);
//                    Vector2D t1 = dot(cl1, nw.tangentVec()) * nw.tangentVec();
//                    Vector2D t2 = dot(cl2, nw.tangentVec()) * nw.tangentVec();

//                    //- m1 is more hydrophobic
//                    GhostCellImmersedBoundaryObject::ZeroGradientStencil m1(cell, *ibObj, cl1);
//                    GhostCellImmersedBoundaryObject::ZeroGradientStencil m2(cell, *ibObj, cl2);
//                    Vector2D grad = BilinearInterpolator(grid_, ibObj->nearestIntersect(cell.centroid())).grad(gamma);

//                    if (m2.bpValue(gamma) < m1.bpValue(gamma))
//                        std::swap(m1, m2);

//                    if (dot(grad, t2) < dot(grad, t1))
//                        std::swap(m1, m2);

//                    if (theta > M_PI_2)
//                    {
//                        m1.init();
//                        eqn.add(m1.cell(), m1.cells(), m1.coeffs());
//                    }
//                    else
//                    {
//                        m2.init();
//                        eqn.add(m2.cell(), m2.cells(), m2.coeffs());
//                    }
//                }

//                for (const Cell &cell: ibObj->solidCells())
//                {
//                    eqn.add(cell, cell, 1.);
//                    eqn.setSource(cell, -0.);
//                }
//            }

//                break;

//            case ImmersedBoundaryObject::HIGH_ORDER:
//            default:
//                throw Exception("SurfaceTensionForce", "contactLineBcs",
//                                "no contact line bcs for immersed boundary type.");
//        }
//    }

//    return eqn;
//}
