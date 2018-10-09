#include "Math/TrilinosAmesosSparseMatrixSolver.h"

#include "DirectForcingImmersedBoundary.h"
#include "DirectForcingImmersedBoundaryLeastSquaresQuadraticStencil.h"

DirectForcingImmersedBoundary::DirectForcingImmersedBoundary(const Input &input,
                                                             const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                             const std::shared_ptr<CellGroup> &domainCells)
    :
      ImmersedBoundary(input, grid, domainCells)
{
    
}

void DirectForcingImmersedBoundary::updateCells()
{
    localIbCells_.clear();
    localSolidCells_.clear();

    cellStatus_->fill(FLUID_CELLS);

    for(auto &ibObj: ibObjs_)
    {
        ibObj->clear();

        for(const Cell &c: ibObj->cellsWithin(*domainCells_))
            if(localSolidCells_.add(c))
            {
                ibObj->addSolidCell(c);
                (*cellStatus_)(c) = SOLID_CELLS;
            }
    }

    cellStatus_->sendMessages();

    for(const Cell &c: *domainCells_)
    {
        if((*cellStatus_)(c) == SOLID_CELLS)
            continue;

        for(const CellLink &nb: c.neighbours())
        {
            if((*cellStatus_)(nb.cell()) == SOLID_CELLS)
            {
                if(localIbCells_.add(c))
                {
                    ibObj(nb.cell())->addIbCell(c);
                    (*cellStatus_)(c) = IB_CELLS;
                }
                break;
            }
        }
    }

    cellStatus_->sendMessages();
}

FiniteVolumeEquation<Vector2D> DirectForcingImmersedBoundary::computeForcingTerm(const VectorFiniteVolumeField &u,
                                                                                 Scalar timeStep,
                                                                                 VectorFiniteVolumeField &fib) const
{
    FiniteVolumeEquation<Vector2D> eqn(fib, 12);

    for(const Cell& cell: u.cells())
    {
        if(localIbCells_.isInSet(cell))
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, *this);
            auto beta = st.interpolationCoeffs(cell.centroid());

            eqn.add(cell, cell, -timeStep);

            int i = 0;
            for(const Cell *cellPtr: st.cells())
            {
                eqn.add(cell, *cellPtr, beta(0, i) * timeStep);
                eqn.addSource(cell, beta(0, i++) * u(*cellPtr));
            }

            for(const auto &cmpt: st.compatPts())
                eqn.addSource(cell, beta(0, i++) * cmpt.velocity());

            for(const Face *facePtr: st.faces())
                switch(u.boundaryType(*facePtr))
                {
                case VectorFiniteVolumeField::FIXED:
                    eqn.addSource(cell, beta(0, i++) * u(*facePtr));
                    break;
                case VectorFiniteVolumeField::SYMMETRY:
                {
                    Vector2D n = facePtr->norm().unitVec();
                    Vector2D t = n.tangentVec();
                    Tensor2D tmp = outer(t, t);
                    eqn.add(cell, cell, beta(0, i) * timeStep * tmp);
                    eqn.addSource(cell, beta(0, i++) * u(cell));
                    break;
                }
                default:
                    throw Exception("DirectForcingImmersedBoundary", "velocityBcs", "grid boundary type not recognized.");
                }

            eqn.addSource(cell, -u(cell));
        }
        else if (localSolidCells_.isInSet(cell))
        {
            eqn.set(cell, cell, -1.);
            eqn.setSource(cell, (ibObj(cell)->velocity(cell.centroid()) - u(cell)) / timeStep);
        }
        else
        {
            eqn.set(cell, cell, -1.);
            eqn.setSource(cell, 0.);
        }
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> DirectForcingImmersedBoundary::computeForcingTerm(const ScalarFiniteVolumeField &rho,
                                                                                 const VectorFiniteVolumeField &u,
                                                                                 Scalar timeStep,
                                                                                 VectorFiniteVolumeField &fib) const
{
    FiniteVolumeEquation<Vector2D> eqn(fib);

    for(const Cell& cell: u.cells())
    {
        if(localIbCells_.isInSet(cell))
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, *this);

            if(st.nReconstructionPoints() < 6)
            {
                throw Exception("DirectForcingImmersedBoundary",
                                "computeForcingTerm",
                                "not enough cells to perform velocity interpolation. Cell id = "
                                + std::to_string(cell.globalId()) + ", proc = " + std::to_string(grid_->comm().rank()));
            }

            Matrix A(st.nReconstructionPoints(), 6);

            for(int i = 0; i < st.cells().size(); ++i)
            {
                Point2D x = st.cells()[i]->centroid();
                A.setRow(i, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
            }

            for(int i = 0; i < st.compatPts().size(); ++i)
            {
                Point2D x = st.compatPts()[i].pt();
                A.setRow(i + st.cells().size(), {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
            }

            Point2D x = cell.centroid();

            Matrix beta = Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * pseudoInverse(A);

            eqn.add(cell, cell, -timeStep / rho(cell));

            for(int i = 0; i < st.cells().size(); ++i)
            {
                eqn.add(cell, *st.cells()[i], beta(0, i) * timeStep / rho(*st.cells()[i]));
                eqn.addSource(cell, beta(0, i) * u(*st.cells()[i]));
            }

            for(int i = 0; i < st.compatPts().size(); ++i)
                eqn.addSource(cell, beta(0, i + st.cells().size()) * st.compatPts()[i].velocity());

            eqn.addSource(cell, -u(cell));
        }
        else if (localSolidCells_.isInSet(cell))
        {
            eqn.set(cell, cell, -timeStep / rho(cell));
            eqn.setSource(cell, ibObj(cell)->velocity(cell.centroid()) - u(cell));
        }
        else
        {
            eqn.set(cell, cell, -1.);
            eqn.setSource(cell, 0.);
        }
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> DirectForcingImmersedBoundary::velocityBcs(VectorFiniteVolumeField &u,
                                                                          const VectorFiniteVolumeField &uTilde,
                                                                          Scalar timeStep) const
{
    FiniteVolumeEquation<Vector2D> eqn(u, 10);

    for(const Cell& cell: u.cells())
    {
        if(localIbCells_.isInSet(cell))
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, *this);
            auto beta = st.interpolationCoeffs(cell.centroid());

            int i = 0;
            for(const Cell *cellPtr: st.cells())
                eqn.add(cell, *cellPtr, beta(0, i++) * cell.volume() / timeStep);

            for(const auto &cmpt: st.compatPts())
                eqn.addSource(cell, beta(0, i++) * cmpt.velocity() * cell.volume() / timeStep);

            for(const Face *facePtr: st.faces())
                switch(u.boundaryType(*facePtr))
                {
                case VectorFiniteVolumeField::FIXED:
                    eqn.addSource(cell, beta(0, i++) * u(*facePtr) * cell.volume() / timeStep);
                    break;
                case VectorFiniteVolumeField::SYMMETRY:
                {
                    Vector2D n = facePtr->norm().unitVec();
                    Vector2D t = n.tangentVec();
                    Tensor2D tmp = outer(t, t);
                    eqn.add(cell, cell, beta(0, i++) * tmp * cell.volume() / timeStep);
                    break;
                }
                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                    eqn.add(cell, cell, beta(0, i++) * cell.volume() / timeStep);
                    break;
                default:
                    throw Exception("DirectForcingImmersedBoundary", "velocityBcs", "grid boundary type not recognized.");
                }

            eqn.addSource(cell, -uTilde(cell) * cell.volume() / timeStep);
        }
        else if (localSolidCells_.isInSet(cell))
            eqn.addSource(cell, (ibObj(cell)->velocity(cell.centroid()) - uTilde(cell)) * cell.volume() / timeStep);
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> DirectForcingImmersedBoundary::polarVelocityBcs(VectorFiniteVolumeField &u,
                                                                               const VectorFiniteVolumeField &uTilde,
                                                                               Scalar timeStep) const
{
    FiniteVolumeEquation<Vector2D> eqn(u, 10);

    for(const Cell& cell: u.cells())
    {
        if(localIbCells_.isInSet(cell))
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, *this);
            auto beta = st.interpolationCoeffs(cell.centroid());
            Scalar vol = cell.polarVolume();
            int i = 0;
            for(const Cell *cellPtr: st.cells())
                eqn.add(cell, *cellPtr, beta(0, i++) * vol / timeStep);

            for(const auto &cmpt: st.compatPts())
                eqn.addSource(cell, beta(0, i++) * cmpt.velocity() * vol / timeStep);

            for(const Face *facePtr: st.faces())
                switch(u.boundaryType(*facePtr))
                {
                case VectorFiniteVolumeField::FIXED:
                    eqn.addSource(cell, beta(0, i++) * u(*facePtr) * vol / timeStep);
                    break;
                case VectorFiniteVolumeField::SYMMETRY:
                {
                    Vector2D n = facePtr->norm().unitVec();
                    Vector2D t = n.tangentVec();
                    Tensor2D tmp = outer(t, t);
                    eqn.add(cell, cell, beta(0, i++) * tmp * vol / timeStep);
                    break;
                }
                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                    eqn.add(cell, cell, beta(0, i++) * vol / timeStep);
                    break;
                default:
                    throw Exception("DirectForcingImmersedBoundary", "polarVelocityBcs", "grid boundary type not recognized.");
                }

            eqn.addSource(cell, -uTilde(cell) * vol / timeStep);
        }
        else if (localSolidCells_.isInSet(cell))
            eqn.addSource(cell, (ibObj(cell)->velocity(cell.centroid()) - uTilde(cell)) * cell.polarVolume() / timeStep);
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> DirectForcingImmersedBoundary::velocityBcs(const ScalarFiniteVolumeField &rho,
                                                                          VectorFiniteVolumeField &u,
                                                                          const VectorFiniteVolumeField &uTilde,
                                                                          Scalar timeStep) const
{
    FiniteVolumeEquation<Vector2D> eqn(u, 10);

    for(const Cell& cell: u.cells())
    {
        if(localIbCells_.isInSet(cell))
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, *this);
            auto beta = st.interpolationCoeffs(cell.centroid());

            int i = 0;
            for(const Cell *cellPtr: st.cells())
                eqn.add(cell, *cellPtr, rho(cell) * beta(0, i++) * cell.volume() / timeStep);

            for(const auto &cmpt: st.compatPts())
                eqn.addSource(cell, rho(cell) * beta(0, i++) * cmpt.velocity() * cell.volume() / timeStep);

            for(const Face *facePtr: st.faces())
                switch(u.boundaryType(*facePtr))
                {
                case VectorFiniteVolumeField::FIXED:
                    eqn.addSource(cell, rho(cell) * beta(0, i++) * u(*facePtr) * cell.volume() / timeStep);
                    break;
                case VectorFiniteVolumeField::SYMMETRY:
                {
                    Vector2D n = facePtr->norm().unitVec();
                    Vector2D t = n.tangentVec();
                    Tensor2D tmp = outer(t, t);
                    eqn.add(cell, cell, rho(cell) * beta(0, i++) * tmp * cell.volume() / timeStep);
                    break;
                }
                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                    eqn.add(cell, cell, rho(cell) * beta(0, i++) * cell.volume() / timeStep);
                    break;
                default:
                    throw Exception("DirectForcingImmersedBoundary", "polarVelocityBcs", "grid boundary type not recognized.");
                }

            eqn.addSource(cell, -rho(cell) * uTilde(cell) * cell.volume() / timeStep);
        }
        else if (localSolidCells_.isInSet(cell))
            eqn.addSource(cell, rho(cell) * (ibObj(cell)->velocity(cell.centroid()) - uTilde(cell)) * cell.volume() / timeStep);
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> DirectForcingImmersedBoundary::polarVelocityBcs(const ScalarFiniteVolumeField &rho,
                                                                               VectorFiniteVolumeField &u,
                                                                               const VectorFiniteVolumeField &uTilde,
                                                                               Scalar timeStep) const
{
    FiniteVolumeEquation<Vector2D> eqn(u, 10);

    for(const Cell& cell: u.cells())
    {
        if(localIbCells_.isInSet(cell))
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, *this);
            auto beta = st.interpolationCoeffs(cell.centroid());
            Scalar vol = cell.polarVolume();
            int i = 0;
            for(const Cell *cellPtr: st.cells())
                eqn.add(cell, *cellPtr, rho(cell) * beta(0, i++) * vol / timeStep);

            for(const auto &cmpt: st.compatPts())
                eqn.addSource(cell, rho(cell) * beta(0, i++) * cmpt.velocity() * vol / timeStep);

            for(const Face *facePtr: st.faces())
                switch(u.boundaryType(*facePtr))
                {
                case VectorFiniteVolumeField::FIXED:
                    eqn.addSource(cell, rho(cell) * beta(0, i++) * u(*facePtr) * vol / timeStep);
                    break;
                case VectorFiniteVolumeField::SYMMETRY:
                {
                    Vector2D n = facePtr->norm().unitVec();
                    Vector2D t = n.tangentVec();
                    Tensor2D tmp = outer(t, t);
                    eqn.add(cell, cell, rho(cell) * beta(0, i++) * tmp * vol / timeStep);
                    break;
                }
                case VectorFiniteVolumeField::NORMAL_GRADIENT:
                    eqn.add(cell, cell, rho(cell) * beta(0, i++) * vol / timeStep);
                    break;
                default:
                    throw Exception("DirectForcingImmersedBoundary", "polarVelocityBcs", "grid boundary type not recognized.");
                }

            eqn.addSource(cell, -rho(cell) * uTilde(cell) * vol / timeStep);
        }
        else if (localSolidCells_.isInSet(cell))
            eqn.addSource(cell, rho(cell) * (ibObj(cell)->velocity(cell.centroid()) - uTilde(cell)) * cell.polarVolume() / timeStep);
    }

    return eqn;
}

void DirectForcingImmersedBoundary::applyHydrodynamicForce(Scalar rho,
                                                           Scalar mu,
                                                           const VectorFiniteVolumeField &u,
                                                           const ScalarFiniteVolumeField &p,
                                                           const Vector2D &g)
{
    struct Stress
    {
        Point2D pt;
        Scalar p;
        Tensor2D tau;
    };

    CrsEquation eqn;
    eqn.setSparseSolver(std::make_shared<TrilinosAmesosSparseMatrixSolver>(grid_->comm(), Tpetra::DynamicProfile));


    for(auto &ibObj: ibObjs_)
    {
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

        eqn.clear();
        Index row = 0;
        for(const Cell &cell: ibObj->ibCells())
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, *this);

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
            tau = mu * (tau + tau.transpose());

            stresses.push_back(Stress{xb, 0., tau});

            auto divTau = mu * Vector2D(4 * coeffs(0, 0) + 2 * coeffs(1, 0) + coeffs(2, 1),
                                        2 * coeffs(0, 1) + 4 * coeffs(1, 1) + coeffs(2, 0));

            eqn.addRows(st.nReconstructionPoints(), 12);

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

            eqn.setRhs(row++, rho * dot(ibObj->acceleration(xb), n) - dot(divTau, n));
        }

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
            stress.p = a * xb.x * xb.x + b * xb.y * xb.y + c * xb.x * xb.y + d * xb.x + e * xb.y + f + rho * dot(xb, g);
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
                }
                else
                {
                    fPressure += (pA + pB) / 2. * (ptA - ptB).normalVec();
                    fShear += dot((tauA + tauB) / 2., (ptB - ptA).normalVec());
                }
            }

            force = fPressure + fShear + ibObj->rho * g * ibObj->shape().area();

            //            std::cout << "Pressure force = " << fPressure << "\n"
            //                      << "Shear force = " << fShear << "\n"
            //                      << "Weight = " << ibObj->rho * g * ibObj->shape().area() << "\n"
            //                      << "Net force = " << force << "\n";
        }

        ibObj->applyForce(grid_->comm().broadcast(grid_->comm().mainProcNo(), force));
    }
}

void DirectForcingImmersedBoundary::applyHydrodynamicForce(const ScalarFiniteVolumeField &rho,
                                                           const ScalarFiniteVolumeField &mu,
                                                           const VectorFiniteVolumeField &u,
                                                           const ScalarFiniteVolumeField &p,
                                                           const Vector2D &g)
{
    struct Stress
    {
        Point2D pt;
        Scalar p;
        Tensor2D tau;
    };

    for(auto &ibObj: ibObjs_)
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
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(cell, *this);

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
            tau = mu(cell) * (tau + tau.transpose());

            stresses.push_back(Stress{xb, 0., tau});

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

void DirectForcingImmersedBoundary::applyHydrodynamicForce(Scalar rho, const VectorFiniteVolumeField &fib, const Vector2D &g)
{
    for(const auto &ibObj: ibObjs_)
    {
        Vector2D f = Vector2D(0., 0.);

        for(const Cell &c: ibObj->cells())
            f -= rho * fib(c) * c.volume();

        f = grid_->comm().sum(f);
        f += (ibObj->rho - rho) * g * ibObj->shape().area();
        ibObj->applyForce(f);
    }
}

void DirectForcingImmersedBoundary::applyHydrodynamicForce(const VectorFiniteVolumeField &fib)
{
    applyHydrodynamicForce(1., fib);
}

//void DirectForcingImmersedBoundary::rce(const ScalarFiniteVolumeField &rho,
//                                                       const ScalarFiniteVolumeField &mu,
//                                                       const VectorFiniteVolumeField &u,
//                                                       const ScalarFiniteVolumeField &p,
//                                                       const ScalarFiniteVolumeField &gamma,
//                                                       const SurfaceTensionForce &ft,
//                                                       const Vector2D &g)
//{
//    Equation eqn;

//    auto nLocalIbCells = grid_->comm().allGather(ibCells_.size());
//    auto colOwnershipRange = std::make_pair(6 * std::accumulate(nLocalIbCells.begin(), nLocalIbCells.begin() + grid_->comm().rank(), 0),
//                                            6 * std::accumulate(nLocalIbCells.begin(), nLocalIbCells.begin() + grid_->comm().rank() + 1, 0));

//    auto colMap = std::vector<Index>(grid_->cells().size(), -1);

//    Index ibCellId = 0;
//    for(const Cell& cell: ibCells_)
//        colMap[cell.id()] = colOwnershipRange.first + 6 * ibCellId++;

//    grid_->sendMessages(colMap);

//    auto ibCells = grid_->globalCellGroup(ibCells_);
//    auto solidCells = grid_->globalCellGroup(solidCells_);

//    std::vector<const Cell*> stCells;
//    std::vector<std::pair<const Cell*, Point2D>> compatPts;

//    Index row = 0;

//    for(const Cell &cell: ibCells_)
//    {
//        stCells.clear();
//        compatPts.clear();

//        for(const CellLink &nb: cell.cellLinks())
//            if(!solidCells.isInGroup(nb.cell()))
//            {
//                stCells.push_back(&nb.cell());

//                if(ibCells.isInGroup(nb.cell()))
//                    compatPts.push_back(std::make_pair(&nb.cell(), nearestIntersect(nb.cell().centroid())));
//            }

//        eqn.setRank(eqn.rank() + stCells.size() + compatPts.size() + 1);

//        Index col = colMap[cell.id()];

//        for(const auto &cellPtr: stCells)
//        {
//            Point2D x = cellPtr->centroid();

//            eqn.setCoeffs(row,
//            {col, col + 1, col + 2, col + 3, col + 4, col + 5},
//            {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});

//            eqn.setRhs(row++, -(p(*cellPtr) + rho(*cellPtr) * dot(g, cellPtr->centroid())));
//        }

//        for(const auto &compatPt: compatPts)
//        {
//            const Cell &cell = *compatPt.first;
//            Point2D x = compatPt.second;

//            eqn.setCoeffs(row,
//            {col, col + 1, col + 2, col + 3, col + 4, col + 5},
//            {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});

//            Index col2 = colMap[cell.id()];

//            eqn.setCoeffs(row++,
//            {col2, col2 + 1, col2 + 2, col2 + 3, col2 + 4, col2 + 5},
//            {-x.x * x.x, -x.y * x.y, -x.x * x.y, -x.x, -x.y, -1.});
//        }

//        Point2D x = cell.centroid();
//        Vector2D n = nearestEdgeNormal(x).unitVec();

//        eqn.setCoeffs(row,
//        {col, col + 1, col + 2, col + 3, col + 4, col + 5},
//        {2. * x.x * n.x, 2. * x.y * n.y, x.y * n.x + x.x * n.y, n.x, n.y, 0.});

//        eqn.setRhs(row++, rho(cell) * (dot(acceleration(x) - g, n)));
//    }

//    eqn.setSparseSolver(std::make_shared<TrilinosAmesosSparseMatrixSolver>(grid_->comm(), Tpetra::DynamicProfile));
//    eqn.setRank(eqn.rank(), 6 * ibCells_.size());
//    eqn.solveLeastSquares();

//    std::vector<std::tuple<Point2D, Scalar, Tensor2D>> stresses;

//    for(int i = 0; i < 6 * ibCells_.size(); i += 6)
//    {
//        Scalar a = eqn.x(i);
//        Scalar b = eqn.x(i + 1);
//        Scalar c = eqn.x(i + 2);
//        Scalar d = eqn.x(i + 3);
//        Scalar e = eqn.x(i + 4);
//        Scalar f = eqn.x(i + 5);

//        Point2D x = nearestIntersect(ibCells_[i / 6].centroid());

//        Scalar pb = a * x.x * x.x + b * x.y * x.y + c * x.x * x.y + d * x.x + e * x.y + f;

//        std::vector<std::pair<Point2D, Vector2D>> pts;
//        pts.reserve(8);

//        const Cell &cell = ibCells_[i / 6];

//        for(const CellLink &nb: cell.cellLinks())
//            if(!solidCells.isInGroup(nb.cell()))
//            {
//                pts.push_back(std::make_pair(nb.cell().centroid(), u(nb.cell())));

//                if(ibCells.isInGroup(nb.cell()))
//                {
//                    auto bp = nearestIntersect(nb.cell().centroid());
//                    pts.push_back(std::make_pair(bp, velocity(bp)));
//                }
//            }

//        auto bp = nearestIntersect(cell.centroid());
//        pts.push_back(std::make_pair(bp, velocity(bp)));

//        Matrix A(pts.size(), 6), rhs(pts.size(), 2);

//        for(int i = 0; i < pts.size(); ++i)
//        {
//            Point2D x = pts[i].first;
//            Point2D u = pts[i].second;

//            A.setRow(i, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
//            rhs.setRow(i, {u.x, u.y});
//        }

//        auto coeffs = solve(A, rhs);
//        auto derivs = Matrix(2, 6, {
//                                 2. * bp.x, 0., bp.y, 1., 0., 0.,
//                                 0., 2. * bp.y, bp.x, 0., 1., 0.
//                             }) * coeffs;

//        //- Then tensor is tranposed here
//        auto tau = Tensor2D(derivs(0, 0), derivs(1, 0), derivs(0, 1), derivs(1, 1));

//        stresses.push_back(std::make_tuple(x, pb, mu(cell) * (tau + tau.transpose())));
//    }

//    stresses = grid_->comm().gatherv(grid_->comm().mainProcNo(), stresses);

//    if(grid_->comm().isMainProc())
//    {

//        std::sort(stresses.begin(), stresses.end(), [this](std::tuple<Point2D, Scalar, Tensor2D> &lhs,
//                  std::tuple<Point2D, Scalar, Tensor2D> &rhs)
//        {
//            return (std::get<0>(lhs) - shape_->centroid()).angle()
//                    < (std::get<0>(rhs) - shape_->centroid()).angle();
//        });

//        Vector2D fShear(0., 0.), fPressure(0., 0.);

//        for(int i = 0; i < stresses.size(); ++i)
//        {
//            auto qpa = stresses[i];
//            auto qpb = stresses[(i + 1) % stresses.size()];

//            auto ptA = std::get<0>(qpa);
//            auto ptB = std::get<0>(qpb);
//            auto pA = std::get<1>(qpa);
//            auto pB = std::get<1>(qpb);
//            auto tauA = std::get<2>(qpa);
//            auto tauB = std::get<2>(qpb);

//            fPressure += (pA + pB) / 2. * (ptA - ptB).normalVec();
//            fShear += dot((tauA + tauB) / 2., (ptB - ptA).normalVec());
//        }

//        force_ = fPressure + fShear + this->rho * g * shape_->area();

//        std::cout << "Pressure force = " << fPressure << "\n"
//                  << "Shear force = " << fShear << "\n"
//                  << "Weight = " << this->rho * g * shape_->area() << "\n"
//                  << "Net force = " << force_ << "\n";
//    }

//    force_ = grid_->comm().broadcast(grid_->comm().mainProcNo(), force_);
//}
