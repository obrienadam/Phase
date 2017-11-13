#include "HighOrderImmersedBoundaryObject.h"
#include "BlockMatrix.h"
#include "HighOrderStencil.h"

HighOrderImmersedBoundaryObject::HighOrderImmersedBoundaryObject(const std::string &name,
                                                                 Label id,
                                                                 FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryObject(name, id, grid)
{

}

void HighOrderImmersedBoundaryObject::updateCells()
{
    fluid_->add(cells_);
    ibCells_.clear();
    solidCells_.clear();
    cells_.addAll(fluid_->itemsWithin(*shapePtr_));
    solidCells_.addAll(cells_);

    for (const Cell &cell: solidCells_)
        for (const InteriorLink &nb: cell.neighbours())
        {
            if (!solidCells_.isInGroup(nb.cell()) && grid().localActiveCells().isInGroup(nb.cell()))
            {
                cells_.add(nb.cell());
                ibCells_.add(nb.cell());
            }
        }

    //constructDirichletCoeffs();
    //constructNeumannCoeffs();
}

Equation<Scalar> HighOrderImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &phi) const
{
    Equation<Scalar> eqn(phi);
    CellGroup compatCells;
    for(const Cell& cell: solidCells_)
    {
        int count = 0;
        for(const CellLink &nb: cell.neighbours())
            if(!isInIb(nb.cell()))
                count++;
        if(count > 1)
            compatCells.add(cell);
    }

    for(const Cell& cell: solidCells_)
        if(!compatCells.isInGroup(cell))
            eqn.add(cell, cell, 1.);

    for(const Cell& cell: ibCells_)
    {
        for(const CellLink& nb: cell.neighbours())
        {
            if(compatCells.isInGroup(nb.cell()))
            {

            }
        }

        for(const CellLink& dg: cell.diagonals())
        {
            //if(!isInIb(dg.cell()))
            //    addFRow()
        }
    }

    return eqn;
}

Equation<Vector2D> HighOrderImmersedBoundaryObject::velocityBcs(VectorFiniteVolumeField &u) const
{
    Equation<Vector2D> eqn(u);

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, -velocity(cell.centroid()));
    }

    for (const Cell &cell: ibCells_)
    {
        std::vector<Ref<const Cell>> cells;
        std::vector<Point2D> bps;

        for(const CellLink& nb: cell.neighbours())
        {
            if(ibCells_.isInGroup(nb.cell()))
            {
                cells.push_back(std::cref(nb.cell()));
                bps.push_back(nearestIntersect(nb.cell().centroid()));
            }
            else if(fluid_->isInGroup(nb.cell()))
                cells.push_back(std::cref(nb.cell()));
        }

        for(const CellLink& dg: cell.diagonals())
        {
            if(ibCells_.isInGroup(dg.cell()))
            {
                cells.push_back(std::cref(dg.cell()));
                bps.push_back(nearestIntersect(dg.cell().centroid()));
            }
            else if(fluid_->isInGroup(dg.cell()))
                cells.push_back(std::cref(dg.cell()));
        }

        bps.push_back(nearestIntersect(cell.centroid()));

        Matrix A(cells.size() + bps.size(), 6);

        auto addRow = [&A](int i, const Point2D& x)
        {
            A(i, 0) = x.x * x.x;
            A(i, 1) = x.y * x.y;
            A(i, 2) = x.x * x.y;
            A(i, 3) = x.x;
            A(i, 4) = x.y;
            A(i, 5) = 1.;
        };

        int i = 0;
        for(const Cell& cell: cells)
            addRow(i++, cell.centroid());

        for(const Point2D& pt: bps)
            addRow(i++, pt);

        Point2D x = cell.centroid();
        auto c = Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * pseudoInverse(A);

        eqn.add(cell, cell, -1.);
        eqn.add(cell, cells.begin(), cells.end(), c.begin());

        i = 0;
        for(auto it = c.begin() + cells.size(); it != c.end(); ++it)
        {
            eqn.addSource(cell, velocity(bps[i++]) * (*it));
        }

        std::cout << A << "\n\n";
        std::cout << c << "\n";
    }
    std::cout << "FINISHED\n";

    return eqn;
}

Equation<Scalar> HighOrderImmersedBoundaryObject::contactLineBcs(ScalarFiniteVolumeField &gamma, Scalar theta) const
{
    Equation<Scalar> eqn(gamma);

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, 0.);
    }

    for (const Cell &cell: ibCells_)
    {
        StaticMatrix<9, 6> A;

        auto addFixedRow = [&A](int i, const Point2D &x) {
            A(i, 0) = x.x * x.x;
            A(i, 1) = x.y * x.y;
            A(i, 2) = x.x * x.y;
            A(i, 3) = x.x;
            A(i, 4) = x.y;
            A(i, 5) = 1.;
        };

        auto addClRow = [&A](int i, const Point2D &bp, const Vector2D &cl) {
            A(i, 0) = 2. * bp.x * cl.x;
            A(i, 1) = 2. * bp.y * cl.y;
            A(i, 2) = bp.x * cl.y + bp.y * cl.x;
            A(i, 3) = cl.x;
            A(i, 4) = cl.y;
            A(i, 5) = 0.;
        };

        auto nbs = cell.cellLinks();
        std::vector<Ref<const Cell>> stCells;
        std::vector<Point2D> bn;

        int i = 0;
        for (const CellLink &nb: nbs)
        {
            if (isInIb(nb.cell()))
                continue;

            addFixedRow(i++, nb.cell().centroid());
            stCells.push_back(nb.cell());
        }

        for (const CellLink &nb: nbs)
        {
            if (!isInIb(nb.cell()))
                continue;

            Point2D bp = nearestIntersect(nb.cell().centroid());
            Vector2D n = nearestEdgeNormal(nb.cell().centroid());
            addClRow(i++, bp, n);
        }


        Point2D x = cell.centroid();
        addClRow(i, x, nearestEdgeNormal(x));

        auto c = StaticMatrix<1, 6>({x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * pseudoInverse(A);

        eqn.add(cell, cell, -1.);
        eqn.add(cell, stCells.begin(), stCells.end(), c.begin());

        for (auto it = c.begin() + stCells.size(); it != c.end(); ++it)
            eqn.addSource(cell, (*it) * 0);
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

        auto addRow = [&A](int i, const Point2D &x) {
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

        auto addFixedRow = [&A](int i, const Point2D &x) {
            A(i, 0) = x.x * x.x;
            A(i, 1) = x.y * x.y;
            A(i, 2) = x.x * x.y;
            A(i, 3) = x.x;
            A(i, 4) = x.y;
            A(i, 5) = 1.;
        };

        auto addDerivRow = [&A](int i, const Point2D &x, const Vector2D &n) {
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