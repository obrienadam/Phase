#include "Math/EigenSparseMatrixSolver.h"

#include "HighOrderImmersedBoundaryObject.h"

HighOrderImmersedBoundaryObject::HighOrderImmersedBoundaryObject(const std::string &name,
                                                                 const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                                 const std::shared_ptr<CellGroup> &solverCells)
        :
        ImmersedBoundaryObject(name, grid, solverCells)
{

}

void HighOrderImmersedBoundaryObject::updateCells()
{
    clear();

    auto items = solverCells_->itemsWithin(*shape_);
    cells_.add(items.begin(), items.end());
    solidCells_.add(items.begin(), items.end());

    for (const Cell &cell: solidCells_)
        for (const InteriorLink &nb: cell.neighbours())
        {
            if (!isInIb(nb.cell().centroid()))
            {
                cells_.add(nb.cell());
                ibCells_.add(nb.cell());
            }
        }

    //constructDirichletCoeffs();
    //constructNeumannCoeffs();
}

FiniteVolumeEquation<Scalar> HighOrderImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &phi) const
{
    typedef Eigen::Triplet<Scalar> Triplet;
    typedef Eigen::SparseMatrix<Scalar> SparseMatrix;
    typedef Eigen::SparseLU<SparseMatrix> Solver;

    FiniteVolumeEquation<Scalar> eqn(phi);

    for (const Cell &cell: solidCells_)
        eqn.add(cell, cell, 1.);

    //- Local ids of ibCells
    std::unordered_map<int, int> idToIndex;
    idToIndex.reserve(ibCells_.size());
    int id = 0;
    for (const Cell &cell: ibCells_)
        idToIndex[cell.id()] = id++;

    int M = 0;
    std::vector<Triplet> triplets;
    std::vector<Ref<const Cell>> cells;

    //- Cell eqns
    for (const Cell &cell: ibCells_)
    {
        for (const CellLink &nb: cell.cellLinks())
            if (!solidCells_.isInGroup(nb.cell()))
            {
                Point2D x = nb.cell().centroid();
                int j = idToIndex[cell.id()] * 6;

                triplets.push_back(Triplet(M, j++, x.x * x.x));
                triplets.push_back(Triplet(M, j++, x.y * x.y));
                triplets.push_back(Triplet(M, j++, x.x * x.y));
                triplets.push_back(Triplet(M, j++, x.x));
                triplets.push_back(Triplet(M, j++, x.y));
                triplets.push_back(Triplet(M, j++, 1.));
                cells.push_back(std::cref(nb.cell()));
                M++;
            }
    }

    //- Compatibility eqns
    for (const Cell &cell: ibCells_)
    {
        for (const CellLink &nb: cell.cellLinks())
        {
            if (ibCells_.isInGroup(nb.cell()))
            {
                Point2D x = nearestIntersect(nb.cell().centroid());
                int j = idToIndex[cell.id()] * 6;
                triplets.push_back(Triplet(M, j++, x.x * x.x));
                triplets.push_back(Triplet(M, j++, x.y * x.y));
                triplets.push_back(Triplet(M, j++, x.x * x.y));
                triplets.push_back(Triplet(M, j++, x.x));
                triplets.push_back(Triplet(M, j++, x.y));
                triplets.push_back(Triplet(M, j++, 1.));

                j = idToIndex[nb.cell().id()] * 6;
                triplets.push_back(Triplet(M, j++, -x.x * x.x));
                triplets.push_back(Triplet(M, j++, -x.y * x.y));
                triplets.push_back(Triplet(M, j++, -x.x * x.y));
                triplets.push_back(Triplet(M, j++, -x.x));
                triplets.push_back(Triplet(M, j++, -x.y));
                triplets.push_back(Triplet(M, j++, -1.));
                M++;
            }
        }
    }

    //- Boundary eqns
    for (const Cell &cell: ibCells_)
    {
        Point2D bp = nearestIntersect(cell.centroid());
        Vector2D bn = nearestEdgeNormal(bp);

        int j = idToIndex[cell.id()] * 6;
//        triplets.push_back(Triplet(M, j++, bp.x * bp.x));
//        triplets.push_back(Triplet(M, j++, bp.y * bp.y));
//        triplets.push_back(Triplet(M, j++, bp.x * bp.y));
//        triplets.push_back(Triplet(M, j++, bp.x));
//        triplets.push_back(Triplet(M, j++, bp.y));
//        triplets.push_back(Triplet(M, j++, 1.));

        triplets.push_back(Triplet(M, j++, 2. * bp.x * bn.x));
        triplets.push_back(Triplet(M, j++, 2. * bp.y * bn.y));
        triplets.push_back(Triplet(M, j++, bp.x * bn.y + bp.y * bn.x));
        triplets.push_back(Triplet(M, j++, bn.x));
        triplets.push_back(Triplet(M, j++, bn.y));
        triplets.push_back(Triplet(M, j++, 0.));
        M++;
    }

    //- Assemble matrix
    SparseMatrix A(M, 6 * ibCells_.size());
    A.setFromTriplets(triplets.begin(), triplets.end());

    triplets.clear();
    for (const Cell &cell: ibCells_)
    {
        int j = idToIndex[cell.id()];
        int i = 6 * j;
        Point2D x = cell.centroid();

        triplets.push_back(Triplet(i++, j, x.x * x.x));
        triplets.push_back(Triplet(i++, j, x.y * x.y));
        triplets.push_back(Triplet(i++, j, x.x * x.y));
        triplets.push_back(Triplet(i++, j, x.x));
        triplets.push_back(Triplet(i++, j, x.y));
        triplets.push_back(Triplet(i++, j, 1.));
    }

    //- X Matrix
    SparseMatrix X(6 * ibCells_.size(), ibCells_.size());
    X.setFromTriplets(triplets.begin(), triplets.end());
    Solver solver;
    SparseMatrix P = A.transpose() * A;
    solver.compute(P);
    X = (solver.solve(X)).transpose() * A.transpose();

    for (const Cell &cell: ibCells_)
    {
        for (int j = 0; j < cells.size(); ++j)
            eqn.add(cell, cells[j], X.coeff(idToIndex[cell.id()], j));

        eqn.add(cell, cell, -1.);
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> HighOrderImmersedBoundaryObject::velocityBcs(VectorFiniteVolumeField &u) const
{
    FiniteVolumeEquation<Vector2D> eqn(u);

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, -velocity(cell.centroid()));
    }

    for (const Cell &cell: ibCells_)
    {
        std::vector<Ref<const Cell>> cells;
        std::vector<Point2D> bps;

        for (const CellLink &nb: cell.neighbours())
        {
            if (ibCells_.isInGroup(nb.cell()))
            {
                cells.push_back(std::cref(nb.cell()));
                bps.push_back(nearestIntersect(nb.cell().centroid()));
            }
            else if (solverCells_->isInGroup(nb.cell()))
                cells.push_back(std::cref(nb.cell()));
        }

        for (const CellLink &dg: cell.diagonals())
        {
            if (ibCells_.isInGroup(dg.cell()))
            {
                cells.push_back(std::cref(dg.cell()));
                bps.push_back(nearestIntersect(dg.cell().centroid()));
            }
            else if (solverCells_->isInGroup(dg.cell()))
                cells.push_back(std::cref(dg.cell()));
        }

        bps.push_back(nearestIntersect(cell.centroid()));

        Matrix A(cells.size() + bps.size(), 6);

        auto addRow = [&A](int i, const Point2D &a, const Point2D &b)
        {
            Point2D x = b;

            A(i, 0) = x.x * x.x;
            A(i, 1) = x.y * x.y;
            A(i, 2) = x.x * x.y;
            A(i, 3) = x.x;
            A(i, 4) = x.y;
            A(i, 5) = 1.;
        };

        int i = 0;
        for (const Cell &kcell: cells)
            addRow(i++, cell.centroid(), kcell.centroid());

        for (const Point2D &pt: bps)
            addRow(i++, cell.centroid(), pt);

        Point2D x = cell.centroid();
        auto c = Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * pseudoInverse(A);

        eqn.add(cell, cell, -1.);
        eqn.add(cell, cells.begin(), cells.end(), c.begin());

        i = 0;
        for (auto it = c.begin() + cells.size(); it != c.end(); ++it)
            eqn.addSource(cell, velocity(bps[i++]) * (*it));
    }

    return eqn;
}

FiniteVolumeEquation<Scalar> HighOrderImmersedBoundaryObject::contactLineBcs(ScalarFiniteVolumeField &gamma, Scalar theta) const
{
    typedef Eigen::Triplet<Scalar> Triplet;
    typedef Eigen::SparseMatrix<Scalar> SparseMatrix;
    typedef Eigen::SparseLU<SparseMatrix> Solver;

    FiniteVolumeEquation<Scalar> eqn(gamma);

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
    }

    //- Local ids of ibCells
    std::unordered_map<int, int> ids;
    ids.reserve(ibCells_.size());
    int id = 0;
    for (const Cell &cell: ibCells_)
        ids[cell.id()] = id++;

    int M = 0;
    std::vector<Triplet> triplets;
    std::vector<Ref<const Cell>> cells;

    //- Cell eqns
    for (const Cell &cell: ibCells_)
    {
        for (const CellLink &nb: cell.cellLinks())
            if (!solidCells_.isInGroup(nb.cell()))
            {
                Point2D x = nb.cell().centroid();
                int j = ids[cell.id()] * 6;

                triplets.push_back(Triplet(M, j++, x.x * x.x));
                triplets.push_back(Triplet(M, j++, x.y * x.y));
                triplets.push_back(Triplet(M, j++, x.x * x.y));
                triplets.push_back(Triplet(M, j++, x.x));
                triplets.push_back(Triplet(M, j++, x.y));
                triplets.push_back(Triplet(M, j++, 1.));
                cells.push_back(std::cref(nb.cell()));
                M++;
            }
    }

    //- Compatibility eqns
    for (const Cell &cell: ibCells_)
    {
        for (const CellLink &nb: cell.cellLinks())
        {
            if (ibCells_.isInGroup(nb.cell()))
            {
                Point2D x = nearestIntersect(nb.cell().centroid());
                int j = ids[cell.id()] * 6;
                triplets.push_back(Triplet(M, j++, x.x * x.x));
                triplets.push_back(Triplet(M, j++, x.y * x.y));
                triplets.push_back(Triplet(M, j++, x.x * x.y));
                triplets.push_back(Triplet(M, j++, x.x));
                triplets.push_back(Triplet(M, j++, x.y));
                triplets.push_back(Triplet(M, j++, 1.));

                j = ids[nb.cell().id()] * 6;
                triplets.push_back(Triplet(M, j++, -x.x * x.x));
                triplets.push_back(Triplet(M, j++, -x.y * x.y));
                triplets.push_back(Triplet(M, j++, -x.x * x.y));
                triplets.push_back(Triplet(M, j++, -x.x));
                triplets.push_back(Triplet(M, j++, -x.y));
                triplets.push_back(Triplet(M, j++, -1.));
                M++;
            }
        }
    }

    //- Boundary eqns
    for (const Cell &cell: ibCells_)
    {
        Point2D bp = nearestIntersect(cell.centroid());
        Vector2D bn = nearestEdgeNormal(bp);

        int j = ids[cell.id()] * 6;
//        triplets.push_back(Triplet(M, j++, bp.x * bp.x));
//        triplets.push_back(Triplet(M, j++, bp.y * bp.y));
//        triplets.push_back(Triplet(M, j++, bp.x * bp.y));
//        triplets.push_back(Triplet(M, j++, bp.x));
//        triplets.push_back(Triplet(M, j++, bp.y));
//        triplets.push_back(Triplet(M, j++, 1.));

        Vector2D cl = bn.rotate(M_PI_2 - theta);
        triplets.push_back(Triplet(M, j++, 2. * bp.x * cl.x));
        triplets.push_back(Triplet(M, j++, 2. * bp.y * cl.y));
        triplets.push_back(Triplet(M, j++, bp.x * cl.y + bp.y * cl.x));
        triplets.push_back(Triplet(M, j++, cl.x));
        triplets.push_back(Triplet(M, j++, cl.y));
        triplets.push_back(Triplet(M++, j++, 0.));

        j = ids[cell.id()] * 6;
        cl = bn.rotate(theta - M_PI_2);
        triplets.push_back(Triplet(M, j++, 2. * bp.x * cl.x));
        triplets.push_back(Triplet(M, j++, 2. * bp.y * cl.y));
        triplets.push_back(Triplet(M, j++, bp.x * cl.y + bp.y * cl.x));
        triplets.push_back(Triplet(M, j++, cl.x));
        triplets.push_back(Triplet(M, j++, cl.y));
        triplets.push_back(Triplet(M++, j++, 0.));
    }

    //- Assemble matrix
    SparseMatrix A(M, 6 * ibCells_.size());
    A.setFromTriplets(triplets.begin(), triplets.end());

    triplets.clear();
    for (const Cell &cell: ibCells_)
    {
        int j = ids[cell.id()];
        int i = 6 * j;
        Point2D x = cell.centroid();

        triplets.push_back(Triplet(i++, j, x.x * x.x));
        triplets.push_back(Triplet(i++, j, x.y * x.y));
        triplets.push_back(Triplet(i++, j, x.x * x.y));
        triplets.push_back(Triplet(i++, j, x.x));
        triplets.push_back(Triplet(i++, j, x.y));
        triplets.push_back(Triplet(i++, j, 1.));
    }

    //- X Matrix
    SparseMatrix X(6 * ibCells_.size(), ibCells_.size());
    X.setFromTriplets(triplets.begin(), triplets.end());
    Solver solver;
    SparseMatrix P = A.transpose() * A;
    solver.compute(P);
    X = (solver.solve(X)).transpose() * A.transpose();

    for (const Cell &cell: ibCells_)
    {
        for (int j = 0; j < cells.size(); ++j)
            eqn.add(cell, cells[j], X.coeff(ids[cell.id()], j));

        eqn.add(cell, cell, -1.);
    }

    return eqn;
}

void HighOrderImmersedBoundaryObject::computeForce(Scalar rho, Scalar mu,
                                                   const VectorFiniteVolumeField &u,
                                                   const ScalarFiniteVolumeField &p,
                                                   const Vector2D &g)
{
    int i = 0;
    for (const Cell &cell: ibCells_)
    {
        StaticMatrix<9, 6> A;

//        auto nbs = cell.cellLinks([this](const CellLink& link){
//            return isInIb(link.cell());
//        });
//
//        nbs = cell.cellLinks()
    }
}

//- private

//void HighOrderImmersedBoundaryObject::constructDirichletCoeffsQuad()
//{
//    Ad_.clear();
//    stDCells_.clear();
//    bd_.clear();
//    bps_.clear();
//
//    for (const Cell &cell: ibCells_)
//    {
//        std::vector<Ref<const Cell>> stCells;
//        std::vector<Point2D> bps;
//        Matrix A(6, 6);
//
//        auto addRow = [&A](int i, const Point2D &x) {
//            A(i, 0) = x.x * x.x;
//            A(i, 1) = x.y * x.y;
//            A(i, 2) = x.x * x.y;
//            A(i, 3) = x.x;
//            A(i, 4) = x.y;
//            A(i, 5) = 1.;
//        };
//
//        auto links = cell.cellLinks();
//        Vector2D bn = nearestEdgeNormal(cell.centroid());
//
//        std::sort(links.begin(), links.end(), [&bn, &cell](const CellLink &l, const CellLink &r) {
//            return dot(l.cell().centroid() - cell.centroid(), bn) < dot(r.cell().centroid() - cell.centroid(), bn);
//        });
//
//        int i = 0;
//        for (auto it = links.begin(); it != links.begin() + 5; ++it)
//        {
//            if (!fluid_->isInGroup(it->get().cell()))
//                continue;
//
//            Vector2D x = it->get().cell().centroid();
//            addRow(i++, x);
//            stCells.push_back(it->get().cell());
//        }
//
//        for (auto it = links.begin(); it != links.begin() + 5; ++it)
//        {
//            if (fluid_->isInGroup(it->get().cell()))
//                continue;
//
//            Vector2D bp = nearestIntersect(it->get().cell().centroid());
//            addRow(i++, bp);
//            bps.push_back(bn);
//        }
//
//        Vector2D bp = nearestIntersect(cell.centroid());
//        addRow(i, bp);
//        bps.push_back(bp);
//
//        Vector2D x = cell.centroid();
//        Ad_.push_back(pseudoInverse(A));
//        bd_.push_back(Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * Ad_.back());
//        stDCells_.push_back(stCells);
//        bps_.push_back(bps);
//    }
//}

void HighOrderImmersedBoundaryObject::constructDirichletCoeffs()
{
    Ad_.clear();
    stDCells_.clear();
    bd_.clear();
    bps_.clear();

    for (const Cell &cell: ibCells_)
    {
        StaticMatrix<9, 6> A;
        std::vector<Ref<const Cell>> stCells;
        std::vector<Point2D> bps;

        auto addRow = [&A](int i, const Point2D &x)
        {
            A(i, 0) = x.x * x.x;
            A(i, 1) = x.y * x.y;
            A(i, 2) = x.x * x.y;
            A(i, 3) = x.x;
            A(i, 4) = x.y;
            A(i, 5) = 1.;
        };

        int i = 0;
        for (const CellLink &nb: cell.cellLinks())
        {
            if (isInIb(nb.cell()))
                continue;

            stCells.push_back(std::cref(nb.cell()));
            addRow(i++, nb.cell().centroid());
        }

        for (const CellLink &nb: cell.cellLinks())
        {
            if (!isInIb(nb.cell()))
                continue;

            bps.push_back(nearestIntersect(nb.cell().centroid()));
            addRow(i++, bps.back());
        }

        bps.push_back(nearestIntersect(cell.centroid()));
        addRow(i, bps.back());

        Point2D x = cell.centroid();

        Ad_.push_back(pseudoInverse(A));
        bd_.push_back(StaticMatrix<1, 6>({x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * Ad_.back());
        stDCells_.push_back(stCells);
        bps_.push_back(bps);
    }
}

//void HighOrderImmersedBoundaryObject::constructNeumannCoeffs()
//{
//    An_.clear();
//    bn_.clear();
//    stNCells_.clear();
//    bns_.clear();
//
//    for (const Cell &cell: ibCells_)
//    {
//        Matrix A(6, 6);
//        std::vector<Ref<const Cell>> stCells;
//        std::vector<Vector2D> bns;
//
//        auto addFixedRow = [&A](int i, const Point2D &x) {
//            A(i, 0) = x.x * x.x;
//            A(i, 1) = x.y * x.y;
//            A(i, 2) = x.x * x.y;
//            A(i, 3) = x.x;
//            A(i, 4) = x.y;
//            A(i, 5) = 1.;
//        };
//
//        auto addDerivRow = [&A](int i, const Point2D &x, const Vector2D &n) {
//            A(i, 0) = 2. * x.x * n.x;
//            A(i, 1) = 2. * x.y * n.y;
//            A(i, 2) = x.y * n.x + x.x * n.y;
//            A(i, 3) = n.x;
//            A(i, 4) = n.y;
//            A(i, 5) = 0.;
//        };
//
//        auto links = cell.cellLinks();
//        Vector2D bn = nearestEdgeNormal(cell.centroid());
//
//        std::sort(links.begin(), links.end(), [&bn, &cell](const CellLink &l, const CellLink &r) {
//            return dot(l.cell().centroid() - cell.centroid(), bn) < dot(r.cell().centroid() - cell.centroid(), bn);
//        });
//
//        int i = 0;
//        for (auto it = links.begin(); it != links.begin() + 5; ++it)
//        {
//            if (isInIb(it->get().cell()))
//                continue;
//
//            Vector2D x = it->get().cell().centroid();
//            addFixedRow(i++, x);
//            stCells.push_back(it->get().cell());
//        }
//
//        for (auto it = links.begin(); it != links.begin() + 5; ++it)
//        {
//            if (!isInIb(it->get().cell()))
//                continue;
//
//            Vector2D bp = nearestIntersect(it->get().cell().centroid());
//            Vector2D bn = nearestEdgeNormal(it->get().cell().centroid());
//            addDerivRow(i++, bp, bn);
//            bns.push_back(bn);
//        }
//
//        Point2D bp = nearestIntersect(cell.centroid());
//        bn = nearestEdgeNormal(cell.centroid());
//        addDerivRow(i, bp, bn);
//        bns.push_back(bn);
//
//        Vector2D x = cell.centroid();
//
//        An_.push_back(inverse(A));
//        bn_.push_back(Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * An_.back());
//        stNCells_.push_back(stCells);
//        bns_.push_back(bns);
//    }
//}

void HighOrderImmersedBoundaryObject::constructNeumannCoeffs()
{
    An_.clear();
    bn_.clear();
    stNCells_.clear();
    bns_.clear();

    for (const Cell &cell: ibCells_)
    {
        StaticMatrix<9, 6> A;
        std::vector<Ref<const Cell>> stCells;
        std::vector<Vector2D> bns;

        auto addFixedRow = [&A](int i, const Point2D &x)
        {
            A(i, 0) = x.x * x.x;
            A(i, 1) = x.y * x.y;
            A(i, 2) = x.x * x.y;
            A(i, 3) = x.x;
            A(i, 4) = x.y;
            A(i, 5) = 1.;
        };

        auto addDerivRow = [&A](int i, const Point2D &x, const Vector2D &n)
        {
            A(i, 0) = 2. * x.x * n.x;
            A(i, 1) = 2. * x.y * n.y;
            A(i, 2) = x.y * n.x + x.x * n.y;
            A(i, 3) = n.x;
            A(i, 4) = n.y;
            A(i, 5) = 0.;
        };

        auto links = cell.cellLinks();

        int i = 0;
        for (const CellLink &nb: links)
        {
            if (isInIb(nb.cell()))
                continue;

            addFixedRow(i++, nb.cell().centroid());
            stCells.push_back(nb.cell());
        }

        for (const CellLink &nb: links)
        {
            if (!isInIb(nb.cell()))
                continue;

            Vector2D bp = nearestIntersect(nb.cell().centroid());
            Vector2D bn = nearestEdgeNormal(nb.cell().centroid());
            addDerivRow(i++, bp, bn);
            bns.push_back(bn);
        }

        Point2D bp = nearestIntersect(cell.centroid());
        Vector2D bn = nearestEdgeNormal(cell.centroid());
        addDerivRow(i, bp, bn);
        bns.push_back(bn);

        Vector2D x = cell.centroid();

        An_.push_back(pseudoInverse(A));
        bn_.push_back(StaticMatrix<1, 6>({x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * An_.back());
        stNCells_.push_back(stCells);
        bns_.push_back(bns);
    }
}