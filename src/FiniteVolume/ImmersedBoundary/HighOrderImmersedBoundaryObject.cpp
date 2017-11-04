#include "HighOrderImmersedBoundaryObject.h"

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

    constructDirichletCoeffs();
}

Equation<Scalar> HighOrderImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &phi) const
{
    Equation<Scalar> eqn(phi);

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, 0.);
    }

    int i = 0;
    for (const Cell &cell: ibCells_)
    {
        const auto &c = bn_[i];
        const auto &stCells = stNCells_[i];

        eqn.add(cell, cell, -1.);
        eqn.add(cell, stCells.begin(), stCells.end(), c.begin());

        for (int j = c.n() - stCells.size(); j < c.n(); ++j)
            eqn.addSource(cell, c(0, j) * 0.); // derivatives

        ++i;
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

    int i = 0;
    for (const Cell &cell: ibCells_)
    {
        const auto &c = bd_[i];
        const auto &stCells = stDCells_[i];

//        std::cout << stCells.size() << "\n"
//                  << c.size() << "\n"
//                  << bps_[i].size() << "\n\n";

        eqn.add(cell, cell, -1.);
        eqn.add(cell, stCells.begin(), stCells.end(), c.begin());

        int j = 0;
        for(auto it = c.begin() + stCells.size(); it != c.end(); ++it)
            eqn.addSource(cell, (*it) * velocity(bps_[i][j++]));

        ++i;
    }

    return eqn;
}

//- private

void HighOrderImmersedBoundaryObject::constructDirichletCoeffs()
{
    Ad_.clear();
    stDCells_.clear();
    bd_.clear();
    bps_.clear();

    for (const Cell &cell: ibCells_)
    {
        Matrix A(9, 6);
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

        Ad_.push_back(inverse(transpose(A) * A) * transpose(A));
        bd_.push_back(Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * Ad_.back());
        stDCells_.push_back(stCells);
        bps_.push_back(bps);
    }
}

void HighOrderImmersedBoundaryObject::constructNeumannCoeffs()
{
    An_.clear();
    bn_.clear();
    stNCells_.clear();
    bns_.clear();

    for (const Cell &cell: ibCells_)
    {
        Matrix A(6, 6);
        std::vector<Ref<const Cell>> stCells;
        std::vector<Vector2D> bns;

        auto addFixedRow = [&A](int i, const Point2D &x) {};

        auto addDerivRow = [&A](int i, const Vector2D &n) {};

        Vector2D x = cell.centroid();

        An_.push_back(inverse(A));
        bn_.push_back(Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * An_.back());
        stNCells_.push_back(stCells);
        bns_.push_back(bns);
    }
}