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
    FiniteVolumeEquation<Scalar> eqn(gamma);

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

void CelesteImmersedBoundary::applyForce(const ScalarFiniteVolumeField &rho,
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
        Scalar p, gamma;
        Tensor2D tau;
    };

    for(auto &ibObj: ib)
    {
        Equation eqn;

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

            //- Compute the surface properties from the contact line
            CelesteImmersedBoundary::ContactLineStencil clst(*ibObj, xb, theta(*ibObj), gamma);

            auto coeffs = solve(A, rhs);
            auto derivs = Matrix(2, 6, {
                                     2. * xb.x, 0., xb.y, 1., 0., 0.,
                                     0., 2. * xb.y, xb.x, 0., 1., 0.
                                 }) * coeffs;

            //- The tensor is tranposed here
            auto tau = Tensor2D(derivs(0, 0), derivs(1, 0), derivs(0, 1), derivs(1, 1));
            tau = mu(cell) * (tau + tau.transpose());

            stresses.push_back(Stress{xb, 0., clst.gamma(), tau});

            auto divTau = mu(cell) * Vector2D(4 * coeffs(0, 0) + 2 * coeffs(1, 0) + coeffs(2, 1),
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

            eqn.setRhs(row++, rho(cell) * dot(ibObj->acceleration(xb), n) - dot(divTau, n));
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

            const Cell &cell = ibObj->ibCells()[i / 6];

            Stress &stress = stresses[i / 6];
            Point2D xb = stress.pt;
            stress.p = a * xb.x * xb.x + b * xb.y * xb.y + c * xb.x * xb.y + d * xb.x + e * xb.y + f + rho(cell) * dot(xb, g);
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

            Vector2D fShear(0., 0.), fPressure(0., 0.);

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

                fPressure += (pA + pB) / 2. * (ptA - ptB).normalVec();
                fShear += dot((tauA + tauB) / 2., (ptB - ptA).normalVec());
            }

            force = fPressure + fShear + ibObj->rho * g * ibObj->shape().area();

            std::cout << "Pressure force = " << fPressure << "\n"
                      << "Shear force = " << fShear << "\n"
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

void CelesteImmersedBoundary::computeCurvature()
{
    kappa_->fill(0.);

    auto &n = *n_;
    auto &kappa = *kappa_;

    for (const Cell &cell: kappa.cells())
    {
        bool isInterface = n(cell).magSqr() > 0. && !ib_.lock()->ibObj(cell.centroid());

        if(isInterface)
            for(const CellLink &nb: cell.cellLinks())
                if(n(nb.cell()).magSqr() == 0.)
                {
                    isInterface = false;
                    break;
                }

        if (isInterface)
            kappa(cell) = kappaStencils_[cell.id()].kappa(n);
    }

    kappa.grid()->sendMessages(kappa);

    for (const Face &face: kappa.grid()->interiorFaces())
    {
        //- According to Afkhami 2007

        if(kappa(face.lCell()) != 0. && kappa(face.rCell()) != 0.)
        {
            Scalar g = face.volumeWeight();
            kappa(face) = g * kappa(face.lCell()) + (1. - g) * kappa(face.rCell());
        }
        else if (kappa(face.lCell()) != 0.)
            kappa(face) = kappa(face.lCell());
        else if (kappa(face.rCell()) != 0.)
            kappa(face) = kappa(face.rCell());
        else
            kappa(face) = 0.;
    }

    for(const Face &face: kappa.grid()->boundaryFaces())
        if(n(face.lCell()).magSqr() != 0.)
            kappa(face) = kappa(face.lCell());
}
