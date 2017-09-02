#include "Celeste.h"
#include "Algorithm.h"
#include "GhostCellImmersedBoundaryObject.h"

Celeste::Celeste(const Input &input,
                 const ImmersedBoundary &ib,
                 const ScalarFiniteVolumeField &gamma,
                 const ScalarFiniteVolumeField &rho,
                 const ScalarFiniteVolumeField &mu,
                 const VectorFiniteVolumeField &u,
                 const ScalarGradient &gradGamma)
        :
        SurfaceTensionForce(input, ib, gamma, rho, mu, u, gradGamma)
{
    constructMatrices();
}

void Celeste::computeFaces()
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    auto &ft = *this;
    auto &kappa = *kappa_;

    for (const Face &face: gamma_.grid().faces())
        ft(face) = sigma_ * kappa(face) * gradGamma_(face);
}

void Celeste::computeFaces(const ImmersedBoundary &ib)
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature(ib);

    auto &ft = *this;
    auto &kappa = *kappa_;

    for (const Face &face: gamma_.grid().faces())
        ft(face) = sigma_ * kappa(face) * gradGamma_(face);
}

void Celeste::compute()
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    auto &ft = *this;
    auto &kappa = *kappa_;

    ft.fill(Vector2D(0., 0.));
    for (const Cell &cell: gamma_.grid().cellZone("fluid"))
        ft(cell) = sigma_ * kappa(cell) * gradGamma_(cell);
}

void Celeste::compute(const ImmersedBoundary &ib)
{
    computeGradGammaTilde(ib);
    computeInterfaceNormals();
    computeCurvature(ib);

    auto &ft = *this;
    const auto &kappa = *kappa_;

    ft.fill(Vector2D(0, 0));

    for (const Cell &cell: grid_->cellZone("fluid"))
        ft(cell) = sigma_ * kappa(cell) * gradGamma_(cell);

    ft.interpolateFaces();
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

    for (const Cell &cell: grid_->localActiveCells())
        gradGammaTildeMatrices_[cell.id()] = leastSquaresMatrix(cell, true);
}

void Celeste::constructKappaMatrices()
{
    kappaMatrices_.resize(grid_->cells().size());

    for (const Cell &cell: grid_->cellZone("fluid"))
        kappaMatrices_[cell.id()] = leastSquaresMatrix(cell, false);
}

void Celeste::computeGradGammaTilde()
{
    smoothGammaField();

    auto &gammaTilde = *gammaTilde_;
    auto &gradGammaTilde = *gradGammaTilde_;

    gradGammaTilde_->fill(Vector2D(0., 0.));
    Matrix b(8, 1);

    for (const Cell &cell: gradGamma_.grid().cellZone("fluid"))
    {
        Size stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        b.resize(stencilSize, 1);

        int i = 0;
        for (const InteriorLink &nb: cell.neighbours())
            b(i++, 0) = (gammaTilde(nb.cell()) - gammaTilde(cell)) / nb.rCellVec().magSqr();

        for (const DiagonalCellLink &dg: cell.diagonals())
            b(i++, 0) = (gammaTilde(dg.cell()) - gammaTilde(cell)) / dg.rCellVec().magSqr();

        for (const BoundaryLink &bd: cell.boundaries())
            b(i++, 0) = (gammaTilde(bd.face()) - gammaTilde(cell)) / bd.rFaceVec().magSqr();

        b = gradGammaTildeMatrices_[cell.id()] * b;
        gradGammaTilde(cell) = Vector2D(b(0, 0), b(1, 0));
    }
}

void Celeste::computeGradGammaTilde(const ImmersedBoundary &ib)
{
    smoothGammaField();

    auto &gammaTilde = *gammaTilde_;

    for (const auto ibObj: ib.ibObjPtrs())
    {
        auto gcIbObj = std::static_pointer_cast<GhostCellImmersedBoundaryObject>(ibObj);
        Scalar theta = getTheta(*gcIbObj);

        for (const GhostCellStencil &st: gcIbObj->stencils())
        {
            gammaTilde(st.cell()) = 0.;

            //- Produce characteristic lines d1 and d2
            Vector2D nb = -gcIbObj->nearestEdgeNormal(st.boundaryPoint());
            Line2D d1(st.cell().centroid(), nb.rotate(M_PI_2 - theta).normalVec());
            Line2D d2(st.cell().centroid(), nb.rotate(theta - M_PI_2).normalVec());

            //- Track whether the line has been constructed
            bool isConstructed[2] = {false, false};
            std::vector<Ref<const Cell>> dCells[2];
            Point2D ip[2];

            const CellZone &fluid = grid_->cellZone("fluid");

            //- Look for stencil candidates
            int nNodes = 1;
            while (!isConstructed[0] || !isConstructed[1])
            {
                //- Find candidate nodes
                for (const Node &node: grid_->findNearestNodes(st.boundaryPoint(), nNodes++))
                {
                    auto cells = node.cells();

                    if (fluid.isInGroup(cells.begin(), cells.end()))
                    {
                        Polygon bilinearStencil = Polygon({
                                cells[0].get().centroid(),
                                cells[1].get().centroid(),
                                cells[2].get().centroid(),
                                cells[3].get().centroid()
                        }).convexHull();

                        if (!isConstructed[0])
                        {
                            ip[0] = d1.nearestPoint(bilinearStencil.centroid());
                            if (bilinearStencil.isInside(ip[0]))
                            {
                                dCells[0] = cells;
                                isConstructed[0] = true;
                            }
                        }

                        if (!isConstructed[1])
                        {
                            ip[1] = d2.nearestPoint(bilinearStencil.centroid());
                            if (bilinearStencil.isInside(ip[1]))
                            {
                                dCells[1] = cells;
                                isConstructed[1] = true;
                            }
                        }
                    }
                }
            }

            //- Stencils now constructed
            for (int i = 0; i < 2; ++i)
            {
                Point2D x[] = {
                        dCells[i][0].get().centroid(),
                        dCells[i][1].get().centroid(),
                        dCells[i][2].get().centroid(),
                        dCells[i][3].get().centroid()
                };

                Matrix A(4, 4, {
                        x[0].x * x[0].y, x[0].x, x[0].y, 1.,
                        x[1].x * x[1].y, x[1].x, x[1].y, 1.,
                        x[2].x * x[2].y, x[2].x, x[2].y, 1.,
                        x[3].x * x[3].y, x[3].x, x[3].y, 1.
                });

                Matrix b(4, 1, {
                        gammaTilde(dCells[i][0]),
                        gammaTilde(dCells[i][1]),
                        gammaTilde(dCells[i][2]),
                        gammaTilde(dCells[i][3])
                });

                gammaTilde(st.cell()) += (Matrix(1, 4, {ip[i].x * ip[i].y, ip[i].x, ip[i].y, 1.}) * A.solve(b))(0, 0);
                gammaTilde(st.cell()) = clamp(gammaTilde(st.cell()), 0., 1.);
            }
        }
    }

    auto &gradGammaTilde = *gradGammaTilde_;
    gradGammaTilde.fill(Vector2D(0., 0.));
    Matrix b(8, 1);

    for (const Cell &cell: gradGamma_.grid().cellZone("fluid"))
    {
        Size stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        b.resize(stencilSize, 1);

        int i = 0;
        for (const InteriorLink &nb: cell.neighbours())
            b(i++, 0) = (gammaTilde(nb.cell()) - gammaTilde(cell)) / nb.rCellVec().magSqr();

        for (const DiagonalCellLink &dg: cell.diagonals())
            b(i++, 0) = (gammaTilde(dg.cell()) - gammaTilde(cell)) / dg.rCellVec().magSqr();

        for (const BoundaryLink &bd: cell.boundaries())
            b(i++, 0) = (gammaTilde(bd.face()) - gammaTilde(cell)) / bd.rFaceVec().magSqr();

        b = gradGammaTildeMatrices_[cell.id()] * b;
        gradGammaTilde(cell) = Vector2D(b(0, 0), b(1, 0));
    }
}

void Celeste::computeCurvature()
{
    Matrix bx(8, 1), by(8, 1);
    auto &n = *n_;
    auto &kappa = *kappa_;

    kappa.fill(0.);
    for (const Cell &cell: grid_->cellZone("fluid"))
    {
        Size stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        bx.resize(stencilSize, 1);
        by.resize(stencilSize, 1);

        if (n(cell) == Vector2D(0., 0.))
            continue;

        bool compute = true;

        int i = 0;
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D dn = n(nb.cell()) - n(cell);
            bx(i, 0) = dn.x;
            by(i++, 0) = dn.y;

            if (n(nb.cell()) == Vector2D(0., 0.))
                compute = false;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            Vector2D dn = n(dg.cell()) - n(cell);
            bx(i, 0) = dn.x;
            by(i++, 0) = dn.y;

            if (n(dg.cell()) == Vector2D(0., 0.))
                compute = false;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D dn = n(bd.face()) - n(cell);
            bx(i, 0) = dn.x;
            by(i++, 0) = dn.y;

            if (n(bd.face()) == Vector2D(0., 0.))
                compute = false;
        }

        if (compute)
        {
            bx = kappaMatrices_[cell.id()] * bx;
            by = kappaMatrices_[cell.id()] * by;

            kappa(cell) = bx(0, 0) + by(1, 0);
        }
    }

    grid_->sendMessages(kappa);
    kappa.interpolateFaces();
}

void Celeste::computeCurvature(const ImmersedBoundary &ib)
{
    Matrix bx(8, 1), by(8, 1);
    kappa_->fill(0.);

    auto &n = *n_;
    auto &kappa = *kappa_;

    if (modifiedStencil_.empty())
        modifiedStencil_.resize(grid_->nCells(), false);

    auto modifyStencil = [&ib](const Cell &cell) -> std::shared_ptr<const ImmersedBoundaryObject> {
        for (const InteriorLink &nb: cell.neighbours())
        {
            auto ibObj = ib.ibObj(nb.cell().centroid());
            if (ibObj)
                return ibObj;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            auto ibObj = ib.ibObj(dg.cell().centroid());
            if (ibObj)
                return ibObj;
        }
        return nullptr;
    };

    for (const Cell &cell: grid_->cellZone("fluid"))
    {
        if (n(cell) == Vector2D(0., 0.))
            continue;

        Size stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        bx.resize(stencilSize, 1);
        by.resize(stencilSize, 1);

        //- Determine if ls stencil needs to be modified
        auto ibObj = modifyStencil(cell);

        if (ibObj)
        {
            kappaMatrices_[cell.id()] = leastSquaresMatrix(*ibObj, cell, false);
            modifiedStencil_[cell.id()] = true;
        }
        else if (modifiedStencil_[cell.id()])
        {
            kappaMatrices_[cell.id()] = leastSquaresMatrix(cell, false);
            modifiedStencil_[cell.id()] = false;
        }

        int i = 0;
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D dn = n(nb.cell()) - n(cell);

            if (ibObj)
                dn = contactLineNormal(cell, nb.cell(), *ibObj) - n(cell);

            bx(i, 0) = dn.x;
            by(i++, 0) = dn.y;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            Vector2D dn = n(dg.cell()) - n(cell);

            if (ibObj)
                dn = contactLineNormal(cell, dg.cell(), *ibObj) - n(cell);

            bx(i, 0) = dn.x;
            by(i++, 0) = dn.y;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D dn = n(bd.face()) - n(cell);

            bx(i, 0) = dn.x;
            by(i++, 0) = dn.y;
        }

        bx = kappaMatrices_[cell.id()] * bx;
        by = kappaMatrices_[cell.id()] * by;

        kappa(cell) = bx(0, 0) + by(1, 0);
    }

    grid_->sendMessages(kappa);
    kappa.interpolateFaces();
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

    return inverse(transpose(A) * A) * transpose(A);
}

Matrix Celeste::leastSquaresMatrix(const ImmersedBoundaryObject &ibObj, const Cell &cell, bool weighted)
{
    Matrix A(
            cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size(),
            5
    );

    auto addRow = [&ibObj, &A, weighted](int row, const Point2D &ptA, const Point2D &ptB) {
        Vector2D r = ibObj.intersectionLine(LineSegment2D(ptA, ptB)).ptB() - ptA;

        //- Soften r
        r = r.magSqr() < (ptB - ptA).magSqr() / 100. ? (ptB - ptA) / 10. : r;

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

    return inverse(transpose(A) * A) * transpose(A);
}

