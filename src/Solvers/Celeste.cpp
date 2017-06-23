#include "Celeste.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"
#include "LineSegment2D.h"

Celeste::Celeste(const Input &input,
                 const ScalarFiniteVolumeField &gamma,
                 const ScalarFiniteVolumeField &rho,
                 const ScalarFiniteVolumeField &mu,
                 const VectorFiniteVolumeField &u,
                 VectorFiniteVolumeField &gradGamma)
        :
        ContinuumSurfaceForce(input, gamma, rho, mu, u, gradGamma)
{
    constructMatrices();
}

void Celeste::compute()
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    auto &ft = *this;
    auto &kappa = *kappa_;

    ft.fill(Vector2D(0., 0.));

    for (const Face &face: gamma_.grid().interiorFaces())
    {
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        gradGamma_(face) = (gamma_(face.rCell()) - gamma_(face.lCell())) * rc / rc.magSqr();
        ft(face) = sigma_ * kappa(face) * gradGamma_(face);
    }

    for (const Cell &cell: gamma_.grid().cellZone("fluid"))
        ft(cell) = sigma_ * kappa(cell) * gradGamma_(cell);

    grid_->sendMessages(ft);
}

void Celeste::compute(const ImmersedBoundary &ib)
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature(ib);

    auto &ft = *this;
    const auto &kappa = *kappa_;

    ft.fill(Vector2D(0, 0));

    for (const Cell &cell: grid_->cellZone("fluid"))
        ft(cell) = sigma_ * kappa(cell) * gradGamma_(cell);

    ft.interpolateFaces([](const Face &face) {
        return face.rCell().volume() / (face.lCell().volume() + face.rCell().volume());
    });

    grid_->sendMessages(ft);
}

void Celeste::constructMatrices()
{
    constructGammaTildeMatrices();
    constructKappaMatrices();
}

//- Protected methods

void Celeste::constructGammaTildeMatrices()
{
    gradGammaTildeMatrices_.resize(grid_->cells().size());

    for(const Cell& cell: grid_->cells())
        gradGammaTildeMatrices_[cell.id()] = leastSquaresMatrix(cell, true);
}

void Celeste::constructKappaMatrices()
{
    kappaMatrices_.resize(grid_->cells().size());

    for(const Cell& cell: grid_->cells())
        kappaMatrices_[cell.id()] = leastSquaresMatrix(cell, false);
}

void Celeste::computeGradGammaTilde()
{
    *gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    grid_->sendMessages(*gammaTilde_);

    //- Used in the reconstruction of gradGammaTilde
    gammaTilde_->setBoundaryFaces();

    auto &gammaTilde = *gammaTilde_;
    auto &gradGammaTilde = *gradGammaTilde_;

    gradGammaTilde_->fill(Vector2D(0., 0.));
    Matrix b(8, 1);

    for (const Cell &cell: gradGamma_.grid().localActiveCells())
    {
        Size stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        b.resize(stencilSize, 1);

        int i = 0;
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar sSqr = nb.rCellVec().magSqr();
            b(i, 0) = (gammaTilde(nb.cell()) - gammaTilde(cell)) / sSqr;
            ++i;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            Scalar sSqr = dg.rCellVec().magSqr();
            b(i, 0) = (gammaTilde(dg.cell()) - gammaTilde(cell)) / sSqr;
            ++i;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar sSqr = bd.rFaceVec().magSqr();
            b(i, 0) = (gammaTilde(bd.face()) - gammaTilde(cell)) / sSqr;
            ++i;
        }

        b = gradGammaTildeMatrices_[cell.id()] * b;
        gradGammaTilde(cell) = Vector2D(b(0, 0), b(1, 0));
    }

    grid_->sendMessages(*gradGammaTilde_);
}

void Celeste::computeCurvature()
{
    Matrix bx(8, 1), by(8, 1);
    kappa_->fill(0.);

    auto &n = *n_;
    auto &gradGammaTilde = *gradGammaTilde_;
    auto &kappa = *kappa_;

    for (const Cell &cell: grid_->cellZone("fluid"))
    {
        if (gradGammaTilde(cell).magSqr() < curvatureCutoffTolerance_)
            continue;

        const size_t stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        int i = 0;
        bx.resize(stencilSize, 1);
        by.resize(stencilSize, 1);

        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D dn = n(nb.cell()) - n(cell);
            Scalar sSqr = pow(nb.rFaceVec().mag() + eps_, 2);
            sSqr = 1;

            bx(i, 0) = dn.x / sSqr;
            by(i, 0) = dn.y / sSqr;

            ++i;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            Vector2D dn = n(dg.cell()) - n(cell);
            Scalar sSqr = pow(dg.rCellVec().mag() + eps_, 2);
            sSqr = 1;

            bx(i, 0) = dn.x / sSqr;
            by(i, 0) = dn.y / sSqr;

            ++i;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D dn = n(bd.face()) - n(cell);
            Scalar sSqr = pow(bd.rFaceVec().mag() + eps_, 2);
            sSqr = 1;

            bx(i, 0) = dn.x / sSqr;
            by(i, 0) = dn.y / sSqr;

            ++i;
        }

        bx = kappaMatrices_[cell.id()] * bx;
        by = kappaMatrices_[cell.id()] * by;

        kappa(cell) = bx(0, 0) + by(1, 0);
    }

    grid_->sendMessages(kappa);
    interpolateCurvatureFaces();
}

void Celeste::computeCurvature(const ImmersedBoundary &ib)
{
    Matrix bx(8, 1), by(8, 1);
    kappa_->fill(0.);

    auto &n = *n_;
    auto &gradGammaTilde = *gradGammaTilde_;
    auto &kappa = *kappa_;

    for (const Cell &cell: grid_->cellZone("fluid"))
    {
        if (gradGammaTilde(cell).magSqr() < curvatureCutoffTolerance_)
            continue;

        const size_t stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        int i = 0;
        bx.resize(stencilSize, 1);
        by.resize(stencilSize, 1);

        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D dn = n(nb.cell()) - n(cell);
            Scalar sSqr = pow(nb.rFaceVec().mag() + eps_, 2);
            sSqr = 1;

            for (const ImmersedBoundaryObject &ibObj: ib.ibObjs())
            {
                if (ibObj.isInIb(nb.cell().centroid()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), nb.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();
                    //sSqr = pow(r.mag() + eps_, 2);

                    dn = computeContactLineNormal(gradGamma_(cell), stencil.second, u_(cell), theta()) -
                         n(cell);

                    n(nb.cell()) = n(cell) + dn;

                    break;
                }
            }

            bx(i, 0) = dn.x / sSqr;
            by(i, 0) = dn.y / sSqr;

            ++i;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            Vector2D dn = n(dg.cell()) - n(cell);
            Scalar sSqr = pow(dg.rCellVec().mag() + eps_, 2);
            sSqr = 1;

            for (const ImmersedBoundaryObject &ibObj: ib.ibObjs())
            {
                if (ibObj.isInIb(dg.cell().centroid()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), dg.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();
                    //sSqr = pow(r.mag() + eps_, 2);
                    std::cout << theta() * 180 / M_PI << std::endl;
                    dn = computeContactLineNormal(gradGammaTilde(cell), stencil.second, u_(cell), theta()) -
                         n(cell);

                    n(dg.cell()) = n(cell) + dn;

                    break;
                }
            }

            bx(i, 0) = dn.x / sSqr;
            by(i, 0) = dn.y / sSqr;

            ++i;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D dn = n(bd.face()) - n(cell);
            Scalar sSqr = pow(bd.rFaceVec().mag() + eps_, 2);
            sSqr = 1;

            bx(i, 0) = dn.x / sSqr;
            by(i, 0) = dn.y / sSqr;

            ++i;
        }

        bx = kappaMatrices_[cell.id()] * bx;
        by = kappaMatrices_[cell.id()] * by;

        kappa(cell) = bx(0, 0) + by(1, 0);
    }

    grid_->sendMessages(kappa);
    interpolateCurvatureFaces();
}

Matrix Celeste::leastSquaresMatrix(const Cell &cell, bool weighted)
{
    Matrix A(
            cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size(),
            5
    );

    auto addRow = [&A, weighted](int row, const Vector2D &r) {
        Scalar w = weighted ? r.magSqr() : 1.;

        A(row, 0) = r.x / w;
        A(row, 1) = r.y / w;
        A(row, 2) = r.x * r.x / (2. * w);
        A(row, 3) = r.y * r.y / (2. * w);
        A(row, 4) = r.x * r.y / w;
    };

    int row = 0;

    for (const InteriorLink &nb: cell.neighbours())
        addRow(row++, nb.rCellVec());

    for (const DiagonalCellLink &dg: cell.diagonals())
        addRow(row++, dg.rCellVec());

    for (const BoundaryLink &bd: cell.boundaries())
        addRow(row++, bd.rFaceVec());

    return A;
}

Matrix Celeste::leastSquaresMatrix(const ImmersedBoundary &ib, const Cell &cell, bool weighted)
{
    Matrix A(
            cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size(),
            5
    );

    auto addRow = [&ib, &A, weighted](int row, const Point2D& ptA, const Point2D &ptB) {
        auto ibObj = ib.ibObj(ptB);
        Vector2D r = ibObj ? ibObj->intersectionLine(ptA, ptB).tan() : ptB - ptA;
        Scalar w = weighted ? r.magSqr() : 1.;

        A(row, 0) = r.x / w;
        A(row, 1) = r.y / w;
        A(row, 2) = r.x * r.x / (2. * w);
        A(row, 3) = r.y * r.y / (2. * w);
        A(row, 4) = r.x * r.y / w;
    };

    int row = 0;

    for (const InteriorLink &nb: cell.neighbours())
        addRow(row++, cell.centroid(), nb.cell().centroid());

    for (const DiagonalCellLink &dg: cell.diagonals())
        addRow(row++, cell.centroid(), dg.cell().centroid());

    for (const BoundaryLink &bd: cell.boundaries())
        addRow(row++, cell.centroid(), bd.face().centroid());

    return A;
}

