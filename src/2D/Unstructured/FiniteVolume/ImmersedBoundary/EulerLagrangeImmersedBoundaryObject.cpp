#include "Math/Matrix.h"
#include "FiniteVolumeGrid2D/StructuredRectilinearGrid.h"

#include "EulerLagrangeImmersedBoundaryObject.h"

EulerLagrangeImmersedBoundaryObject::EulerLagrangeImmersedBoundaryObject(const std::string &name,
                                                                         const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                                         const std::shared_ptr<CellGroup> &solverCells)
        :
        ImmersedBoundaryObject(name, grid, solverCells)
{
    auto eqGrid = std::dynamic_pointer_cast<const StructuredRectilinearGrid>(grid_);

    if (eqGrid && eqGrid->isEquidistant())
    {
        h_ = eqGrid->h();
    }
    else
        throw Exception("EulerLagrangeImmersedBoundaryObject",
                        "EulerLagrangeImmersedBoundaryObject",
                        "must use an equidistant grid.");
}

void EulerLagrangeImmersedBoundaryObject::initCircle(const Point2D &center, Scalar radius)
{
    ImmersedBoundaryObject::initCircle(center, radius);

    Scalar ds = 1.2 * h_;
    initLagrangePoints(shape_->perimeter() / ds + 0.5);
}

void EulerLagrangeImmersedBoundaryObject::updateCells()
{
    updateLagrangePoints();
}

FiniteVolumeEquation<Vector2D> EulerLagrangeImmersedBoundaryObject::velocityBcs(VectorFiniteVolumeField &u) const
{
    FiniteVolumeEquation<Vector2D> eqn(u);

    Matrix D(lagrangePoints_.size(), ibCells_.size());

    for (int i = 0; i < lagrangePoints_.size(); ++i)
        for (int j = 0; j < ibCells_.size(); ++j)
            D(i, j) = kernel(lagrangePoints_[i], ibCells_[j].centroid());

    Matrix Ainv = inverse(h_ * h_ * D * transpose(D));
    Matrix ul(lagrangePoints_.size(), 2);

    for (int i = 0; i < lagrangePoints_.size(); ++i)
    {
        ul(i, 0) = velocity(lagrangePoints_[i]).x;
        ul(i, 1) = velocity(lagrangePoints_[i]).y;
    }

    Matrix rhs = transpose(D) * Ainv * ul;
    Matrix coeffs = transpose(D) * Ainv * (-h_ * h_) * D;

    for (int i = 0; i < ibCells_.size(); ++i)
    {
        for (int j = 0; j < ibCells_.size(); ++j)
            eqn.add(ibCells_[i], ibCells_[j], coeffs(i, j));
        eqn.addSource(ibCells_[i], Vector2D(rhs(i, 0), rhs(i, 1)));
    }

    return eqn;
}

void EulerLagrangeImmersedBoundaryObject::correctVelocity(VectorFiniteVolumeField &u) const
{
    Matrix D(lagrangePoints_.size(), ibCells_.size());

    for (int i = 0; i < lagrangePoints_.size(); ++i)
        for (int j = 0; j < ibCells_.size(); ++j)
            D(i, j) = kernel(lagrangePoints_[i], ibCells_[j].centroid());

    Matrix A = h_ * h_ * D * transpose(D);
    Matrix ul(lagrangePoints_.size(), 2), us(ibCells_.size(), 2);

    for (int i = 0; i < lagrangePoints_.size(); ++i)
    {
        ul(i, 0) = velocity(lagrangePoints_[i]).x;
        ul(i, 1) = velocity(lagrangePoints_[i]).y;
    }

    for (int j = 0; j < ibCells_.size(); ++j)
    {
        us(j, 0) = u(ibCells_[j]).x;
        us(j, 1) = u(ibCells_[j]).y;
    }

    Matrix a = solve(A, ul - h_ * h_ * D * us);
    Matrix du = transpose(D) * a;

    for (int j = 0; j < ibCells_.size(); ++j)
        u(ibCells_[j]) += Vector2D(du(j, 0), du(j, 1));
}

Scalar EulerLagrangeImmersedBoundaryObject::kernel(const Point2D &x, const Point2D &xl) const
{
    auto dr = [](Scalar r)
    {
        r = std::abs(r);

//        if (r <= 1.)
//            return (3. - 2. * r + std::sqrt(1. + 4. * r - 4. * r * r)) / 8.;
//        else if (r <= 2.)
//            return (5. - 2. * r - std::sqrt(-7. + 12. * r - 4. * r * r)) / 8.;
//        else
//            return 0.;

//        if (r < 1.)
//            return 2. / 3. - r * r + r * r * r / 2.;
//        else if (r < 2.)
//            return 4. / 3. - 2. * r + r * r - r * r * r / 6.;
//        else
//            return 0.;

        if (r < 1.)
            return 11. / 20. - std::pow(r, 2) / 2. + std::pow(r, 4) / 4. - std::pow(r, 5) / 12.;
        else if (r < 2.)
            return 17. / 40. + 5. * r / 8. - 7. * std::pow(r, 2) / 4 + 5. * std::pow(r, 3) / 4. - 3 * std::pow(r, 4) / 8
                   + std::pow(r, 5) / 24;
        else if (r < 3)
            return 81. / 40. - 27. * r / 8. + 9. * std::pow(r, 2) / 4. - 3. * std::pow(r, 3) / 4. + std::pow(r, 4) / 8.
                   - std::pow(r, 5) / 120.;
        else
            return 0.;
    };

    Vector2D r = x - xl;
    return dr(r.x / h_) * dr(r.y / h_) / (h_ * h_);
}

//- Private

void EulerLagrangeImmersedBoundaryObject::initLagrangePoints(int nLagrangePoints)
{
    lagrangePoints_.clear();
    lagrangeStencils_.clear();
    ibCells_.clear();

    switch (shape_->type())
    {
        case Shape2D::CIRCLE:
        {
            Scalar dTheta = 2. * M_PI / nLagrangePoints;
            Scalar radius = static_cast<const Circle &>(*shape_).radius();

            for (int i = 0; i < nLagrangePoints; ++i)
            {
                Point2D pt = shape_->centroid() + radius * Vector2D(std::cos(i * dTheta), std::sin(i * dTheta));
                lagrangePoints_.push_back(pt);

                if (solverCells_)
                {
                    lagrangeStencils_.push_back(
                            solverCells_->itemsCoveredBy(Box(
                                    pt - 3. * Vector2D(h_, h_),
                                    pt + 3. * Vector2D(h_, h_)
                            ))
                    );

                    ibCells_.add(lagrangeStencils_.back().begin(), lagrangeStencils_.back().end());
                }
            }
        }
            break;
        default:
            throw Exception("EulerLagrangeImmersedBoundaryObject", "initLagrangePoints", "shape is not supported.");
    }
}