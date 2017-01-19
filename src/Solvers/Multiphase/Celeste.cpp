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

VectorFiniteVolumeField Celeste::compute()
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    VectorFiniteVolumeField ft(gamma_.grid, "ft");

    for(const Face &face: gamma_.grid.faces())
        ft(face) = sigma_*kappa_(face)*gradGamma_(face);

    for(const Cell& cell: gamma_.grid.cellZone("fluid"))
    {
        Scalar sumSfx = 0., sumSfy = 0.;
        ft(cell) = Vector2D(0., 0.);

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D& sf = nb.outwardNorm();

            ft(cell) += Vector2D(ft(nb.face()).x*fabs(sf.x), ft(nb.face()).y*fabs(sf.y))/rho_(nb.face());

            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D& sf = bd.outwardNorm();

            ft(cell) += Vector2D(ft(bd.face()).x*fabs(sf.x), ft(bd.face()).y*fabs(sf.y))/rho_(bd.face());

            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        ft(cell) = rho_(cell)*Vector2D(ft(cell).x/sumSfx, ft(cell).y/sumSfy);
    }

    return ft;
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
    gradGammaTildeStencils_.resize(gradGammaTilde_.size());
    for(const Cell &cell: gammaTilde_.grid.cellZone("fluid"))
    {
        gradGammaTildeStencils_[cell.id()] = gradGammaTilde_.grid.activeCells().kNearestNeighbourSearch(cell.centroid(), 9);

        A.resize(8, 5);
        int i = 0;
        for(const Cell &kCell: gradGammaTildeStencils_[cell.id()])
        {
            if(&cell == &kCell)
                continue;

            const Scalar sSqr = (kCell.centroid() - cell.centroid()).magSqr();
            const Scalar dx = kCell.centroid().x - cell.centroid().x;
            const Scalar dy = kCell.centroid().y - cell.centroid().y;

            A(i, 0) = dx/sSqr;
            A(i, 1) = dy/sSqr;
            A(i, 2) = dx*dx/(2.*sSqr);
            A(i, 3) = dy*dy/(2.*sSqr);
            A(i, 4) = dx*dy/sSqr;

            ++i;
        }

        gradGammaTildeMatrices_[cell.id()] = inverse(transpose(A)*A)*transpose(A);
    }
}

void Celeste::constructKappaMatrices()
{
    Matrix A(8, 5);

    kappaMatrices_.resize(gamma_.size());
    for(const Cell &cell: gamma_.grid.cellZone("fluid"))
    {
        const size_t stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();

        A.resize(stencilSize, 5);

        int i = 0;
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar dx = nb.cell().centroid().x - cell.centroid().x;
            Scalar dy = nb.cell().centroid().y - cell.centroid().y;
            Scalar sSqr = (nb.cell().centroid() - cell.centroid()).magSqr();

            for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
            {
                if(ibObj.cells().isInGroup(nb.cell()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), nb.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();

                    sSqr = std::max(1e-8, r.magSqr());
                    dx = r.x;
                    dy = r.y;

                    break;
                }
            }

            A(i, 0) = dx/sSqr;
            A(i, 1) = dy/sSqr;
            A(i, 2) = dx*dx/(2.*sSqr);
            A(i, 3) = dy*dy/(2.*sSqr);
            A(i, 4) = dx*dy/sSqr;

            ++i;
        }

        for(const DiagonalCellLink &dg: cell.diagonals())
        {
            Scalar dx = dg.cell().centroid().x - cell.centroid().x;
            Scalar dy = dg.cell().centroid().y - cell.centroid().y;
            Scalar sSqr = (dg.cell().centroid() - cell.centroid()).magSqr();

            for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
            {
                if(ibObj.cells().isInGroup(dg.cell()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), dg.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();

                    sSqr = std::max(1e-8, r.magSqr());
                    dx = r.x;
                    dy = r.y;

                    break;
                }
            }

            A(i, 0) = dx/sSqr;
            A(i, 1) = dy/sSqr;
            A(i, 2) = dx*dx/(2.*sSqr);
            A(i, 3) = dy*dy/(2.*sSqr);
            A(i, 4) = dx*dy/sSqr;

            ++i;
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar dx = bd.face().centroid().x - cell.centroid().x;
            Scalar dy = bd.face().centroid().y - cell.centroid().y;
            Scalar sSqr = (bd.face().centroid() - cell.centroid()).magSqr();

            A(i, 0) = dx/sSqr;
            A(i, 1) = dy/sSqr;
            A(i, 2) = dx*dx/(2.*sSqr);
            A(i, 3) = dy*dy/(2.*sSqr);
            A(i, 4) = dx*dy/sSqr;

            ++i;
        }

        kappaMatrices_[cell.id()] = inverse(transpose(A)*A)*transpose(A);
    }
}

void Celeste::computeGradGammaTilde()
{
    gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    gradGammaTilde_.fill(Vector2D(0., 0.));

    Matrix b(8, 1);

    for(const Cell &cell: gradGammaTilde_.grid.cellZone("fluid"))
    {
        const auto &nb = gradGammaTildeStencils_[cell.id()];

        b.resize(8, 1);
        int i = 0;

        for(const Cell &kCell: nb)
        {
            if(&cell == &kCell)
                continue;

            const Scalar sSqr = (kCell.centroid() - cell.centroid()).magSqr();
            b(i, 0) = (gammaTilde_(kCell) - gammaTilde_(cell))/sSqr;

            ++i;
        }

        b = gradGammaTildeMatrices_[cell.id()]*b;
        gradGammaTilde_(cell) = Vector2D(b(0, 0), b(1, 0));
    }
}

void Celeste::computeInterfaceNormals()
{
    for(const Cell &cell: n_.grid.cellZone("fluid"))
        n_(cell) = gradGammaTilde_(cell) == Vector2D(0., 0.) ? Vector2D(0., 0.) : -gradGammaTilde_(cell).unitVec();

    for(const Face &face: n_.grid.interiorFaces()) // Not super important how this is computed, it's just for the cicsam scheme
    {
        const Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        n_(face) = -(gammaTilde_(face.rCell()) - gammaTilde_(face.lCell()))*rc/rc.magSqr();
        n_(face) = n_(face).magSqr() < 1e-15 ? Vector2D(0., 0.) : n_(face).unitVec();
    }

    for(const Face &face: n_.grid.boundaryFaces())
    {
        if(isContactLinePatch(face.patch()))
            n_(face) = computeContactLineNormal(gradGammaTilde_(face.lCell()), face.outwardNorm(face.lCell().centroid()), u_(face.lCell()));
        else
            n_(face) = n_(face.lCell());
    }
}

void Celeste::computeCurvature()
{
    Matrix bx(8, 1), by(8, 1);
    kappa_.fill(0.);

    for(const Cell &cell: kappa_.grid.cellZone("fluid"))
    {
        if(gradGamma_(cell).magSqr() < curvatureCutoffTolerance_)
            continue;

        const size_t stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();
        int i = 0;
        bx.resize(stencilSize, 1);
        by.resize(stencilSize, 1);

        for(const InteriorLink &nb: cell.neighbours())
        {
            Vector2D n = n_(nb.cell()) - n_(cell);
            Scalar sSqr = (nb.cell().centroid() - cell.centroid()).magSqr();

            for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
            {
                if(ibObj.cells().isInGroup(nb.cell()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), nb.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();
                    sSqr = std::max(1e-8, r.magSqr());
                    n = computeContactLineNormal(gradGammaTilde_(cell), stencil.second, u_(cell), ibTheta(ibObj)) - n_(cell);

                    break;
                }
            }

            bx(i, 0) = n.x/sSqr;
            by(i, 0) = n.y/sSqr;

            ++i;
        }

        for(const DiagonalCellLink &dg: cell.diagonals())
        {
            Vector2D n = n_(dg.cell()) - n_(cell);
            Scalar sSqr = (dg.cell().centroid() - cell.centroid()).magSqr();

            for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
            {
                if(ibObj.cells().isInGroup(dg.cell()))
                {
                    auto stencil = ibObj.intersectionStencil(cell.centroid(), dg.cell().centroid());

                    Vector2D r = stencil.first - cell.centroid();
                    sSqr = std::max(1e-8, r.magSqr());

                    n = computeContactLineNormal(gradGammaTilde_(cell), stencil.second, u_(cell), ibTheta(ibObj)) - n_(cell);
                    break;
                }
            }

            bx(i, 0) = n.x/sSqr;
            by(i, 0) = n.y/sSqr;

            ++i;
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D n = n_(bd.face()) - n_(cell);
            Scalar sSqr = (bd.face().centroid() - cell.centroid()).magSqr();

            bx(i, 0) = n.x/sSqr; // contact line faces should already be computed
            by(i, 0) = n.y/sSqr;

            ++i;
        }

        bx = kappaMatrices_[cell.id()]*bx;
        by = kappaMatrices_[cell.id()]*by;

        kappa_(cell) = bx(0, 0) + by(1, 0);
    }

    //weightCurvatures();
    interpolateCurvatureFaces();
}


void Celeste::weightCurvatures()
{
    const auto pow8 = [](Scalar x){ return x*x*x*x*x*x*x*x; };
    const Scalar cutoff = 1e-12;

    for(const Cell &cell: w_.grid.cellZone("fluid"))
        w_(cell) = pow8(1. - 2*fabs(0.5 - std::max(0., std::min(1., gammaTilde_(cell)))));

    kappa_.savePreviousIteration();

    for(const Cell &cell: kappa_.grid.cellZone("fluid"))
    {
        Scalar sumKappaW = kappa_.prevIter()(cell)*w_(cell), sumW = w_(cell);

        for(const InteriorLink& nb: cell.neighbours())
        {
            sumKappaW += kappa_.prevIter()(nb.cell())*w_(nb.cell());
            sumW += w_(nb.cell());
        }

        for(const DiagonalCellLink& dg: cell.diagonals())
        {
            sumKappaW += kappa_.prevIter()(dg.cell())*w_(dg.cell());
            sumW += w_(dg.cell());
        }

        if(sumW < cutoff)
            kappa_(cell) = 0.;
        else
            kappa_(cell) = sumKappaW/sumW;
    }

    kappa_.savePreviousIteration();

    for(const Cell &cell: kappa_.grid.cellZone("fluid"))
    {
        const ScalarFiniteVolumeField& kappaPrev = kappa_.prevIter();

        Scalar num = kappaPrev(cell)*w_(cell), den = w_(cell);

        for(const InteriorLink& nb: cell.neighbours())
        {
            if(gradGammaTilde()(nb.cell()).magSqr() < cutoff)
                continue;

            Scalar w = w_(nb.cell())*pow8(dot(n_(cell), nb.rCellVec().unitVec()));

            num += kappaPrev(nb.cell())*w;
            den += w;
        }

        for(const DiagonalCellLink& dg: cell.diagonals())
        {
            if(gradGammaTilde()(dg.cell()).magSqr() < cutoff)
                continue;

            Scalar w = w_(dg.cell())*pow8(dot(n_(cell), dg.rCellVec().unitVec()));

            num += kappaPrev(dg.cell())*w;
            den += w;
        }

        if(den < cutoff)
            kappa_(cell) = 0.;
        else
            kappa_(cell) = num/den;
    }
 }
