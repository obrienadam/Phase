#include "Celeste.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"
#include "LineSegment2D.h"

Celeste::Celeste(const Input &input,
                 Solver &solver)
    :
      ContinuumSurfaceForce(input, solver),
      w_(solver.addScalarField("w"))
{
    constructMatrices();
}

void Celeste::compute(VectorFiniteVolumeField& ft)
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    ft.fill(Vector2D(0., 0.));

    for (const Face &face: gamma_.grid().interiorFaces())
    {
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        gradGamma_(face) = (gamma_(face.rCell()) - gamma_(face.lCell())) * rc / rc.magSqr();

        bool lCellIsInIb = solver_.ib().isIbCell(face.lCell());
        bool rCellIsInIb = solver_.ib().isIbCell(face.rCell());

        if (!(lCellIsInIb || rCellIsInIb))
            ft(face) = sigma_ * kappa_(face) * gradGamma_(face);
        else
            ft(face) = Vector2D(0., 0.);
    }

    for (const Cell &cell: gamma_.grid().cellZone("fluid"))
    {
        Scalar sumSfx = 0., sumSfy = 0.;
        ft(cell) = Vector2D(0., 0.);

        for (const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D &sf = nb.outwardNorm();

            ft(cell) += Vector2D(ft(nb.face()).x * fabs(sf.x), ft(nb.face()).y * fabs(sf.y)) /
                    rho_(nb.face());

            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D &sf = bd.outwardNorm();

            ft(cell) += Vector2D(ft(bd.face()).x * fabs(sf.x), ft(bd.face()).y * fabs(sf.y)) /
                    rho_(bd.face());

            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        ft(cell) = rho_(cell) * Vector2D(ft(cell).x / sumSfx, ft(cell).y / sumSfy);
    }

    solver().grid().sendMessages(solver().comm(), ft);
}

void Celeste::constructMatrices()
{
    constructGammaTildeMatrices();
    constructKappaMatrices();
}

//- Protected methods

void Celeste::constructGammaTildeMatrices()
{
    Matrix A(8, 5);

    gradGammaTildeMatrices_.resize(gradGammaTilde_.size());

    for (const Cell &cell: gradGamma_.grid().cells())
    {
        Size stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        A.resize(stencilSize, 5);

        int i = 0;
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar sSqr = nb.rCellVec().magSqr();
            Scalar dx = nb.rCellVec().x;
            Scalar dy = nb.rCellVec().y;

            A(i, 0) = dx / sSqr;
            A(i, 1) = dy / sSqr;
            A(i, 2) = dx * dx / (2. * sSqr);
            A(i, 3) = dy * dy / (2. * sSqr);
            A(i, 4) = dx * dy / sSqr;

            ++i;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            Scalar sSqr = dg.rCellVec().magSqr();
            Scalar dx = dg.rCellVec().x;
            Scalar dy = dg.rCellVec().y;

            A(i, 0) = dx / sSqr;
            A(i, 1) = dy / sSqr;
            A(i, 2) = dx * dx / (2. * sSqr);
            A(i, 3) = dy * dy / (2. * sSqr);
            A(i, 4) = dx * dy / sSqr;

            ++i;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar sSqr = bd.rFaceVec().magSqr();
            Scalar dx = bd.rFaceVec().x;
            Scalar dy = bd.rFaceVec().y;

            A(i, 0) = dx / sSqr;
            A(i, 1) = dy / sSqr;
            A(i, 2) = dx * dx / (2. * sSqr);
            A(i, 3) = dy * dy / (2. * sSqr);
            A(i, 4) = dx * dy / sSqr;

            ++i;
        }

        gradGammaTildeMatrices_[cell.id()] = inverse(transpose(A) * A) * transpose(A);
    }
}

void Celeste::constructKappaMatrices()
{
    Matrix A(8, 5);

    kappaMatrices_.resize(kappa_.size());
    for (const Cell &cell: kappa_.grid().cellZone("fluid"))
    {
        Size stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        A.resize(stencilSize, 5);

        int i = 0;
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar dx = nb.rCellVec().x;
            Scalar dy = nb.rCellVec().y;
            Scalar sSqr = pow(nb.rCellVec().mag() + eps_, 2);
            sSqr = 1.;

            for (const ImmersedBoundaryObject &ibObj: solver_.ibObjs())
            {
                if (ibObj.isInIb(nb.cell().centroid()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), nb.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();
                    //sSqr = pow(r.mag() + eps_, 2);

                    dx = r.x;
                    dy = r.y;

                    break;
                }
            }

            A(i, 0) = dx / sSqr;
            A(i, 1) = dy / sSqr;
            A(i, 2) = dx * dx / (2. * sSqr);
            A(i, 3) = dy * dy / (2. * sSqr);
            A(i, 4) = dx * dy / sSqr;

            ++i;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            Scalar dx = dg.rCellVec().x;
            Scalar dy = dg.rCellVec().y;
            Scalar sSqr = pow(dg.rCellVec().mag() + eps_, 2);
            sSqr = 1;

            for (const ImmersedBoundaryObject &ibObj: solver_.ibObjs())
            {
                if (ibObj.isInIb(dg.cell().centroid()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), dg.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();
                    //sSqr = pow(r.mag() + eps_, 2);

                    dx = r.x;
                    dy = r.y;

                    break;
                }
            }

            A(i, 0) = dx / sSqr;
            A(i, 1) = dy / sSqr;
            A(i, 2) = dx * dx / (2. * sSqr);
            A(i, 3) = dy * dy / (2. * sSqr);
            A(i, 4) = dx * dy / sSqr;

            ++i;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar dx = bd.rFaceVec().x;
            Scalar dy = bd.rFaceVec().y;
            Scalar sSqr = pow(bd.rFaceVec().mag() + eps_, 2);
            sSqr = 1;

            A(i, 0) = dx / sSqr;
            A(i, 1) = dy / sSqr;
            A(i, 2) = dx * dx / (2. * sSqr);
            A(i, 3) = dy * dy / (2. * sSqr);
            A(i, 4) = dx * dy / sSqr;

            ++i;
        }

        kappaMatrices_[cell.id()] = inverse(transpose(A) * A) * transpose(A);
    }
}

void Celeste::computeGradGammaTilde()
{
    gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    solver().grid().sendMessages(solver().comm(), gammaTilde_);

    //- Used in the reconstruction of gradGammaTilde
    gammaTilde_.setBoundaryFaces();

    gradGammaTilde_.fill(Vector2D(0., 0.));
    Matrix b(8, 1);

    for (const Cell &cell: gradGamma_.grid().localActiveCells())
    {
        Size stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        b.resize(stencilSize, 1);

        int i = 0;
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar sSqr = nb.rCellVec().magSqr();
            b(i, 0) = (gammaTilde_(nb.cell()) - gammaTilde_(cell)) / sSqr;
            ++i;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            Scalar sSqr = dg.rCellVec().magSqr();
            b(i, 0) = (gammaTilde_(dg.cell()) - gammaTilde_(cell)) / sSqr;
            ++i;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar sSqr = bd.rFaceVec().magSqr();
            b(i, 0) = (gammaTilde_(bd.face()) - gammaTilde_(cell)) / sSqr;
            ++i;
        }

        b = gradGammaTildeMatrices_[cell.id()] * b;
        gradGammaTilde_(cell) = Vector2D(b(0, 0), b(1, 0));
    }

    solver().grid().sendMessages(solver().comm(), gradGammaTilde_);
}

void Celeste::computeInterfaceNormals()
{
    //for(const Cell &cell: n_.grid.cells())
    //    n_(cell) = gradGammaTilde_(cell).magSqr() < curvatureCutoffTolerance_ ? Vector2D(0., 0.) : -gradGammaTilde_(cell).unitVec().unitVec();

    for (const Cell &cell: n_.grid().cells())
        n_(cell) = gradGammaTilde_(cell).magSqr() < curvatureCutoffTolerance_ ? Vector2D(0., 0.) : -gradGammaTilde_(
                                                                                    cell).unitVec();

    solver().grid().sendMessages(solver().comm(), n_);


    for(const Patch& patch: n_.grid().patches())
    {
        if(isContactLinePatch(patch))
        {
            for(const Face& face: patch)
            {
                Vector2D sf = face.outwardNorm(face.lCell().centroid());
                n_(face) = computeContactLineNormal(gradGammaTilde_(face.lCell()), sf, u_(face.lCell()), theta());
            }
        }
        else
        {
            for(const Face& face: patch)
                n_(face) = n_(face.lCell());
        }
    }
}

void Celeste::computeCurvature()
{
    Matrix bx(8, 1), by(8, 1);
    kappa_.fill(0.);

    for (const Cell &cell: kappa_.grid().cellZone("fluid"))
    {
        if (gradGammaTilde_(cell).magSqr() < curvatureCutoffTolerance_)
            continue;

        const size_t stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        int i = 0;
        bx.resize(stencilSize, 1);
        by.resize(stencilSize, 1);

        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D dn = n_(nb.cell()) - n_(cell);
            Scalar sSqr = pow(nb.rFaceVec().mag() + eps_, 2);
            sSqr = 1;

            for (const ImmersedBoundaryObject &ibObj: solver_.ibObjs())
            {
                if (ibObj.isInIb(nb.cell().centroid()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), nb.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();
                    //sSqr = pow(r.mag() + eps_, 2);

                    dn = computeContactLineNormal(gradGamma_(cell), stencil.second, u_(cell), ibTheta(ibObj)) -
                            n_(cell);

                    n_(nb.cell()) = n_(cell) + dn;

                    break;
                }
            }

            bx(i, 0) = dn.x / sSqr;
            by(i, 0) = dn.y / sSqr;

            ++i;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            Vector2D dn = n_(dg.cell()) - n_(cell);
            Scalar sSqr = pow(dg.rCellVec().mag() + eps_, 2);
            sSqr = 1;

            for (const ImmersedBoundaryObject &ibObj: solver_.ibObjs())
            {
                if (ibObj.isInIb(dg.cell().centroid()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), dg.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();
                    //sSqr = pow(r.mag() + eps_, 2);

                    dn = computeContactLineNormal(gradGammaTilde_(cell), stencil.second, u_(cell), ibTheta(ibObj)) -
                            n_(cell);

                    n_(dg.cell()) = n_(cell) + dn;

                    break;
                }
            }

            bx(i, 0) = dn.x / sSqr;
            by(i, 0) = dn.y / sSqr;

            ++i;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D dn = n_(bd.face()) - n_(cell);
            Scalar sSqr = pow(bd.rFaceVec().mag() + eps_, 2);
            sSqr = 1;

            bx(i, 0) = dn.x / sSqr;
            by(i, 0) = dn.y / sSqr;

            ++i;
        }

        bx = kappaMatrices_[cell.id()] * bx;
        by = kappaMatrices_[cell.id()] * by;

        kappa_(cell) = bx(0, 0) + by(1, 0);
    }

    //weightCurvatures();

    solver().grid().sendMessages(solver().comm(), kappa_);
    interpolateCurvatureFaces();
}


void Celeste::weightCurvatures()
{
    Scalar cutoff = 1e-10;

    for (const Cell &cell: w_.grid().cellZone("fluid"))
        w_(cell) = pow(0.5 * (cos(M_PI * (2 * gamma_(cell) - 1)) + 1), 2);

    kappa_.savePreviousIteration();

    for (const Cell &cell: kappa_.grid().cellZone("fluid"))
    {
        Scalar sumKappaW = kappa_.savePreviousIteration()(cell) * w_(cell), sumW = w_(cell);

        for (const InteriorLink &nb: cell.neighbours())
        {
            sumKappaW += kappa_.prevIteration()(nb.cell()) * w_(nb.cell());
            sumW += w_(nb.cell());
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            sumKappaW += kappa_.prevIteration()(dg.cell()) * w_(dg.cell());
            sumW += w_(dg.cell());
        }

        if (sumW < cutoff)
            kappa_(cell) = 0.;
        else
            kappa_(cell) = sumKappaW / sumW;
    }

    solver().grid().sendMessages(solver().comm(), kappa_);
    return;

    kappa_.savePreviousIteration();

    for (const Cell &cell: kappa_.grid().cellZone("fluid"))
    {
        const ScalarFiniteVolumeField &kappaPrev = kappa_.prevIteration();

        Scalar num = kappaPrev(cell) * w_(cell), den = w_(cell);

        for (const InteriorLink &nb: cell.neighbours())
        {
            if (gradGammaTilde()(nb.cell()).magSqr() < cutoff)
                continue;

            Scalar w = w_(nb.cell()) * pow(dot(n_(cell), nb.rCellVec().unitVec()), 8);

            num += kappaPrev(nb.cell()) * w;
            den += w;
        }

        for (const DiagonalCellLink &dg: cell.diagonals())
        {
            if (gradGammaTilde()(dg.cell()).magSqr() < cutoff)
                continue;

            Scalar w = w_(dg.cell()) * pow(dot(n_(cell), dg.rCellVec().unitVec()), 8);

            num += kappaPrev(dg.cell()) * w;
            den += w;
        }

        if (den < cutoff)
            kappa_(cell) = 0.;
        else
            kappa_(cell) = num / den;
    }

    solver().grid().sendMessages(solver().comm(), kappa_);
}
