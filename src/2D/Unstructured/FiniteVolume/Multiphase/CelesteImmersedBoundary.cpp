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
                    ContactLineStencil st(cell,
                                          *ibObj,
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
                    ContactLineStencil st(cell,
                                          *ibObj,
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

        Index row = 0;
        for(const Cell &cell: ibObj->ibCells())
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, ib);

            eqn.setRank(eqn.rank() + st.nReconstructionPoints() + 1);

            Index col = cellIdToIndexMap[cell.id()];

            for(const auto &cellPtr: st.cells())
            {
                Point2D x = cellPtr->centroid();

                eqn.setCoeffs(row,
                {col, col + 1, col + 2, col + 3, col + 4, col + 5},
                {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});

                eqn.setRhs(row++, -p(*cellPtr));
            }

            for(const auto &compatPt: st.compatPts())
            {
                if(&cell == &compatPt.cell())
                    continue;

                Point2D x = compatPt.pt();

                eqn.setCoeffs(row,
                {col, col + 1, col + 2, col + 3, col + 4, col + 5},
                {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});

                Index col2 = cellIdToIndexMap[compatPt.cell().id()];

                eqn.setCoeffs(row++,
                {col2, col2 + 1, col2 + 2, col2 + 3, col2 + 4, col2 + 5},
                {-x.x * x.x, -x.y * x.y, -x.x * x.y, -x.x, -x.y, -1.});
            }

            Point2D x = ibObj->nearestIntersect(cell.centroid());
            Vector2D n = ibObj->nearestEdgeUnitNormal(x);

            eqn.setCoeffs(row,
            {col, col + 1, col + 2, col + 3, col + 4, col + 5},
            {2. * x.x * n.x, 2. * x.y * n.y, x.y * n.x + x.x * n.y, n.x, n.y, 0.});

            eqn.setRhs(row++, rho(cell) * dot(ibObj->acceleration(x), n));
        }

        eqn.setSparseSolver(std::make_shared<TrilinosAmesosSparseMatrixSolver>(grid_->comm(), Tpetra::DynamicProfile));
        eqn.setRank(eqn.rank(), 6 * ibObj->ibCells().size());
        eqn.solveLeastSquares();

        std::vector<Stress> stresses;
        stresses.reserve(ibObj->ibCells().size());

        for(int i = 0; i < 6 * ibObj->ibCells().size(); i += 6)
        {
            Scalar a = eqn.x(i);
            Scalar b = eqn.x(i + 1);
            Scalar c = eqn.x(i + 2);
            Scalar d = eqn.x(i + 3);
            Scalar e = eqn.x(i + 4);
            Scalar f = eqn.x(i + 5);

            const Cell &cell = ibObj->ibCells()[i / 6];

            Stress stress;
            stress.pt = ibObj->nearestIntersect(cell.centroid());

            Point2D x = stress.pt;
            stress.p = a * x.x * x.x + b * x.y * x.y + c * x.x * x.y + d * x.x + e * x.y + f + rho(cell) * dot(x, g);

            std::vector<std::pair<Point2D, Vector2D>> pts;
            pts.reserve(8);

            for(const CellLink &nb: cell.cellLinks())
                if(!ib.globalSolidCells().isInSet(nb.cell()))
                {
                    pts.push_back(std::make_pair(nb.cell().centroid(), u(nb.cell())));

                    if(ib.globalIbCells().isInSet(nb.cell()))
                    {
                        auto bp = ibObj->nearestIntersect(nb.cell().centroid());
                        pts.push_back(std::make_pair(bp, ibObj->velocity(bp)));
                    }
                }

            pts.push_back(std::make_pair(stress.pt, ibObj->velocity(stress.pt)));

            Matrix A(pts.size(), 6), rhs(pts.size(), 2);

            for(int i = 0; i < pts.size(); ++i)
            {
                Point2D x = pts[i].first;
                Point2D u = pts[i].second;

                A.setRow(i, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
                rhs.setRow(i, {u.x, u.y});
            }

            auto coeffs = solve(A, rhs);
            auto derivs = Matrix(2, 6, {
                                     2. * x.x, 0., x.y, 1., 0., 0.,
                                     0., 2. * x.y, x.x, 0., 1., 0.
                                 }) * coeffs;

            //- Then tensor is tranposed here
            auto tau = Tensor2D(derivs(0, 0), derivs(1, 0), derivs(0, 1), derivs(1, 1));

            stress.tau = mu(cell) * (tau + tau.transpose());
            stresses.push_back(stress);
        }

        stresses = grid_->comm().gatherv(grid_->comm().mainProcNo(), stresses);

        Vector2D force(0., 0.);

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
                    ContactLineStencil st(cell,
                                          *ibObj,
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
