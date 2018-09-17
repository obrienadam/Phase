#include "Math/TrilinosAmesosSparseMatrixSolver.h"

#include "CelesteImmersedBoundary.h"
#include "Geometry/Intersection.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundaryLeastSquaresQuadraticStencil.h"
#include "FiniteVolume/ImmersedBoundary/SurfaceField.h"

CelesteImmersedBoundary::CelesteImmersedBoundary(const Input &input,
                                                 const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                 const std::shared_ptr<CellGroup> &fluidCells,
                                                 const std::weak_ptr<const ImmersedBoundary> &ib)
    :
      Celeste(input, grid, fluidCells),
      ib_(ib)
{
    for(const auto& ibObj: *ib_.lock())
    {
        ibContactAngles_[ibObj->name()] = input.boundaryInput().get<Scalar>(
                    "ImmersedBoundaries." + ibObj->name() + ".gamma.contactAngle",
                    90.) * M_PI / 180.;
    }
}

Scalar CelesteImmersedBoundary::theta(const ImmersedBoundaryObject &ibObj) const
{
    auto it = ibContactAngles_.find(ibObj.name());
    return it != ibContactAngles_.end() ? it->second : M_PI_2;
}

void CelesteImmersedBoundary::computeContactLineExtension(ScalarFiniteVolumeField &gamma) const
{
    for(const std::shared_ptr<ImmersedBoundaryObject> &ibObj: *ib_.lock())
        for(const Cell& cell: ibObj->cells())
            if(ibObj->isInIb(cell.centroid()))
            {
                Scalar distSqr = (ibObj->nearestIntersect(cell.centroid()) - cell.centroid()).magSqr();

                if(distSqr <= kernelWidth_ * kernelWidth_)
                {
                    ContactLineStencil st(*ibObj,
                                          cell.centroid(),
                                          ibContactAngles_.find(ibObj->name())->second,
                                          gamma);

                    gamma(cell) = st.gamma();
                }
            }
}

FiniteVolumeEquation<Scalar> CelesteImmersedBoundary::contactLineBcs(ScalarFiniteVolumeField &gamma, Scalar timeStep) const
{
    FiniteVolumeEquation<Scalar> eqn(gamma, 7);

    for(const std::shared_ptr<ImmersedBoundaryObject> &ibObj: *ib_.lock())
        for(const Cell& cell: ibObj->cells())
        {
            Scalar theta = ibContactAngles_.find(ibObj->name())->second;

            if(ibObj->isInIb(cell.centroid()))
            {
                Scalar distSqr = (ibObj->nearestIntersect(cell.centroid()) - cell.centroid()).magSqr();

                if(distSqr <= kernelWidth_ * kernelWidth_)
                {
                    ContactLineStencil st(*ibObj,
                                          cell.centroid(),
                                          theta,
                                          gamma);

                    Scalar alpha = st.link().alpha(st.cl()[2]);

                    //- Add a source
                    eqn.add(cell, st.link().self(), alpha * cell.volume() / timeStep);
                    eqn.add(cell, st.link().cell(), (1. - alpha) * cell.volume() / timeStep);
                    eqn.add(cell, cell, -cell.volume() / timeStep);
                }
            }
        }

    return eqn;
}

void CelesteImmersedBoundary::appyFluidForces(const ScalarFiniteVolumeField &rho,
                                              const ScalarFiniteVolumeField &mu,
                                              const VectorFiniteVolumeField &u,
                                              const ScalarFiniteVolumeField &p,
                                              const ScalarFiniteVolumeField &gamma,
                                              const Vector2D &g,
                                              DirectForcingImmersedBoundary &ib) const
{
    struct Stress
    {
        Point2D pt;
        Scalar rho, p, gamma;
        Vector2D tcl;
        Tensor2D tau;
    };

    for(auto &ibObj: ib)
    {
        CrsEquation eqn;

        auto nLocalCells = grid_->comm().allGather(ibObj->ibCells().size());

        auto indexStart = 6 * std::accumulate(
                    nLocalCells.begin(),
                    nLocalCells.begin() + grid_->comm().rank(), 0);

        auto cellIdToIndexMap = std::vector<Index>(grid_->cells().size(), -1);

        Index ibCellId = 0;
        for(const Cell& cell: ibObj->ibCells())
            cellIdToIndexMap[cell.id()] = indexStart + 6 * ibCellId++;

        grid_->sendMessages(cellIdToIndexMap);

        std::vector<Stress> stresses;
        stresses.reserve(ibObj->ibCells().size());

        Index row = 0;
        for(const Cell &cell: ibObj->ibCells())
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, ib);

            //- Compute the stress tensor
            Matrix A(st.nReconstructionPoints(), 6), rhs(st.nReconstructionPoints(), 2);

            int i = 0;
            for(const Cell *cell: st.cells())
            {
                Point2D x = cell->centroid();
                A.setRow(i, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
                rhs.setRow(i++, {u(*cell).x, u(*cell).y});
            }

            for(const auto &compatPt: st.compatPts())
            {
                Point2D x = compatPt.pt();
                A.setRow(i, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
                rhs.setRow(i++, {compatPt.velocity().x, compatPt.velocity().y});
            }

            Point2D xb = ibObj->nearestIntersect(cell.centroid());
            auto coeffs = solve(A, rhs);
            auto derivs = Matrix(2, 6, {
                                     2. * xb.x, 0., xb.y, 1., 0., 0.,
                                     0., 2. * xb.y, xb.x, 0., 1., 0.
                                 }) * coeffs;

            //- The tensor is tranposed here

            auto tau = Tensor2D(derivs(0, 0), derivs(1, 0), derivs(0, 1), derivs(1, 1));


            auto clst = ContactLineStencil(*ibObj,
                                           xb,
                                           theta(*ibObj),
                                           gamma);

            Scalar mub = clst.interpolate(mu);

            tau = mub * (tau + tau.transpose());

            stresses.push_back(Stress{xb, clst.interpolate(rho), 0., clst.gamma(), clst.tcl(), tau});

            auto divTau = mub * Vector2D(4 * coeffs(0, 0) + 2 * coeffs(1, 0) + coeffs(2, 1),
                                         2 * coeffs(0, 1) + 4 * coeffs(1, 1) + coeffs(2, 0));

            eqn.setRank(eqn.rank() + st.nReconstructionPoints());

            Index colStart = cellIdToIndexMap[cell.id()];
            for(const Cell *cell: st.cells())
            {
                Point2D x = cell->centroid();

                eqn.setCoeffs(row,
                {colStart, colStart + 1, colStart + 2, colStart + 3, colStart + 4, colStart + 5},
                {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});

                eqn.setRhs(row++, -p(*cell));
            }

            for(const auto &compatPt: st.compatPts())
            {
                if(&cell == &compatPt.cell())
                    continue;

                Point2D x = compatPt.pt();

                eqn.setCoeffs(row,
                {colStart, colStart + 1, colStart + 2, colStart + 3, colStart + 4, colStart + 5},
                {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});

                Index colStart2 = cellIdToIndexMap[compatPt.cell().id()];

                eqn.setCoeffs(row++,
                {colStart2, colStart2 + 1, colStart2 + 2, colStart2 + 3, colStart2 + 4, colStart2 + 5},
                {-x.x * x.x, -x.y * x.y, -x.x * x.y, -x.x, -x.y, -1.});
            }

            Vector2D n = ibObj->nearestEdgeUnitNormal(xb);

            eqn.setCoeffs(row,
            {colStart, colStart + 1, colStart + 2, colStart + 3, colStart + 4, colStart + 5},
            {2. * xb.x * n.x, 2. * xb.y * n.y, xb.y * n.x + xb.x * n.y, n.x, n.y, 0.});

            eqn.setRhs(row++, stresses.back().rho * dot(ibObj->acceleration(xb), n) - dot(divTau, n));
        }

        eqn.setSparseSolver(std::make_shared<TrilinosAmesosSparseMatrixSolver>(grid_->comm(), Tpetra::DynamicProfile));
        eqn.setRank(eqn.rank(), 6 * ibObj->ibCells().size());
        eqn.solveLeastSquares();

        //- Extract the pressure coeffs and compute
        for(int i = 0; i < 6 * ibObj->ibCells().size(); i += 6)
        {
            Scalar a = eqn.x(i);
            Scalar b = eqn.x(i + 1);
            Scalar c = eqn.x(i + 2);
            Scalar d = eqn.x(i + 3);
            Scalar e = eqn.x(i + 4);
            Scalar f = eqn.x(i + 5);

            Stress &stress = stresses[i / 6];
            Point2D xb = stress.pt;
            Scalar rhob = stress.rho;

            stress.p = a * xb.x * xb.x + b * xb.y * xb.y + c * xb.x * xb.y + d * xb.x + e * xb.y + f + rhob * dot(xb, g);
        }

        stresses = grid_->comm().gatherv(grid_->comm().mainProcNo(), stresses);

        Vector2D force(0., 0.);

        //- Integrate the stresses
        if(grid_->comm().isMainProc())
        {

            std::sort(stresses.begin(), stresses.end(), [&ibObj](const Stress &lhs, const Stress &rhs)
            {
                return (lhs.pt - ibObj->shape().centroid()).angle()
                        < (rhs.pt - ibObj->shape().centroid()).angle();
            });

            Vector2D fShear(0., 0.), fPressure(0., 0.), fCapillary(0., 0.);

            for(int i = 0; i < stresses.size(); ++i)
            {
                auto qpa = stresses[i];
                auto qpb = stresses[(i + 1) % stresses.size()];

                auto ptA = qpa.pt;
                auto ptB = qpb.pt;
                auto pA = qpa.p;
                auto pB = qpb.p;
                auto tauA = qpa.tau;
                auto tauB = qpb.tau;
                auto gA = qpa.gamma;
                auto gB = qpb.gamma;

                if(ibObj->shape().type() == Shape2D::CIRCLE)
                {
                    const Circle &c = static_cast<const Circle&>(ibObj->shape());
                    Scalar r = c.radius();
                    Scalar tA = (ptA - c.centroid()).angle();
                    Scalar tB = (ptB - c.centroid()).angle();

                    while(tB < tA)
                        tB += 2 * M_PI;

                    fPressure += r * Vector2D(
                                pA*tA*sin(tA) - pA*tB*sin(tA) + pA*cos(tA) - pA*cos(tB) - pB*tA*sin(tB) + pB*tB*sin(tB) - pB*cos(tA) + pB*cos(tB),
                                -pA*tA*cos(tA) + pA*tB*cos(tA) + pA*sin(tA) - pA*sin(tB) + pB*tA*cos(tB) - pB*tB*cos(tB) - pB*sin(tA) + pB*sin(tB)
                                ) / (tA - tB);

                    fShear += r * Vector2D(
                                -tA*tauA.xx*sin(tA) + tA*tauA.xy*cos(tA) + tA*tauB.xx*sin(tB) - tA*tauB.xy*cos(tB) + tB*tauA.xx*sin(tA) - tB*tauA.xy*cos(tA) - tB*tauB.xx*sin(tB) + tB*tauB.xy*cos(tB) - tauA.xx*cos(tA) + tauA.xx*cos(tB) - tauA.xy*sin(tA) + tauA.xy*sin(tB) + tauB.xx*cos(tA) - tauB.xx*cos(tB) + tauB.xy*sin(tA) - tauB.xy*sin(tB),
                                -tA*tauA.yx*sin(tA) + tA*tauA.yy*cos(tA) + tA*tauB.yx*sin(tB) - tA*tauB.yy*cos(tB) + tB*tauA.yx*sin(tA) - tB*tauA.yy*cos(tA) - tB*tauB.yx*sin(tB) + tB*tauB.yy*cos(tB) - tauA.yx*cos(tA) + tauA.yx*cos(tB) - tauA.yy*sin(tA) + tauA.yy*sin(tB) + tauB.yx*cos(tA) - tauB.yx*cos(tB) + tauB.yy*sin(tA) - tauB.yy*sin(tB)
                                ) / (tA - tB);

                    // Scalar dgamma = std::abs(gB - gA) / (tB - tA);
                    Scalar theta = this->theta(*ibObj);

                    //- sharp method
                    if(gA < 0.5 != gB <= 0.5)
                    {
                        Scalar alpha = (0.5 - gB) / (gA - gB);
                        Scalar phi = alpha * tA + (1. - alpha) * tB;

                        Vector2D t1 = Vector2D(std::cos(phi - M_PI_2 + theta), std::sin(phi - M_PI_2 + theta));
                        Vector2D t2 = Vector2D(std::cos(phi + M_PI_2 - theta), std::sin(phi + M_PI_2 - theta));

                        fCapillary += sigma_ * (dot(t1, qpa.tcl) > dot(t2, qpa.tcl) ? t1 : t2);
                    }

                    // - diffuse method

                    //                    Vector2D t1 = Vector2D(std::cos(tA - M_PI_2 + theta), std::sin(tA - M_PI_2 + theta));
                    //                    Vector2D t2 = Vector2D(std::cos(tA + M_PI_2 - theta), std::sin(tA + M_PI_2 - theta));

                    //                    if(dot(t1, qpa.tcl) > dot(t2, qpa.tcl))
                    //                    {
                    //                        fCapillary += Vector2D(
                    //                                    dgamma*sigma_*(cos(tA + theta) - cos(tB + theta)),
                    //                                    dgamma*sigma_*(sin(tA + theta) - sin(tB + theta))
                    //                                    );
                    //                    }
                    //                    else
                    //                    {
                    //                        fCapillary += Vector2D(
                    //                                    dgamma*sigma_*(-cos(tA - theta) + cos(tB - theta)),
                    //                                    dgamma*sigma_*(-sin(tA - theta) + sin(tB - theta))
                    //                                    );
                    //                    }
                }
                else
                {
                    fPressure += (pA + pB) / 2. * (ptA - ptB).normalVec();
                    fShear += dot((tauA + tauB) / 2., (ptB - ptA).normalVec());

                    if(gA < 0.5 != gB <= 0.5)
                    {
                        Scalar alpha = (0.5 - gB) / (gA - gB);
                        Vector2D tcl = (alpha * qpa.tcl + (1. - alpha) * qpb.tcl).unitVec();
                        fCapillary += sigma_ * tcl;
                    }
                }
            }

            force = fPressure + fShear + fCapillary + ibObj->rho * g * ibObj->shape().area();

            std::cout << "Pressure force = " << fPressure << "\n"
                      << "Shear force = " << fShear << "\n"
                      << "Capillary force = " << fCapillary << "\n"
                      << "Weight = " << ibObj->rho * g * ibObj->shape().area() << "\n"
                      << "Net force = " << force << "\n";
        }

        ibObj->applyForce(grid_->comm().broadcast(grid_->comm().mainProcNo(), force));
    }
}

void CelesteImmersedBoundary::computeInterfaceNormals()
{
    const VectorFiniteVolumeField &gradGammaTilde = *gradGammaTilde_;
    VectorFiniteVolumeField &n = *n_;

    for (const Cell &cell: n.grid()->cells())
        n(cell) = gradGammaTilde(cell).magSqr() >= eps_ * eps_ ? -gradGammaTilde(cell).unitVec() : Vector2D(0., 0.);

    //- Override the ib cells
    for(const auto &ibObj: *ib_.lock())
        for(const Cell& cell: ibObj->cells())
            if(ibObj->isInIb(cell.centroid()))
            {
                Scalar distSqr = (ibObj->nearestIntersect(cell.centroid()) - cell.centroid()).magSqr();

                if(distSqr <= kernelWidth_ * kernelWidth_)
                {
                    ContactLineStencil st(*ibObj,
                                          cell.centroid(),
                                          ibContactAngles_.find(ibObj->name())->second,
                                          *gammaTilde_);

                    n(cell) = st.ncl();
                }
            }

    //- Boundary faces set from contact line orientation
    for (const FaceGroup &patch: n.grid()->patches())
    {
        Scalar theta = SurfaceTensionForce::theta(patch);

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
