#include "Celeste.h"

Celeste::Celeste(const Input &input,
                 Solver &solver)
    :
      ContinuumSurfaceForce(input, solver),
      wGamma_(solver.addScalarField("wGamma"))
{
    constructMatrices();
}

VectorFiniteVolumeField Celeste::compute()
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    VectorFiniteVolumeField ft(gamma_.grid, "ft");
    VectorFiniteVolumeField gradGamma = grad(gamma_);

    for(const Cell &cell: gamma_.grid.fluidCells())
        ft[cell.id()] = sigma_*kappa_[cell.id()]*gradGamma[cell.id()];

    return ft;
}

void Celeste::constructMatrices()
{
    Matrix A(8, 5);

    gradGammaTildeMatrices_.resize(gradGammaTilde_.size());
    gradGammaTildeStencils_.resize(gradGammaTilde_.size());
    for(const Cell &cell: gammaTilde_.grid.fluidCells())
    {
        gradGammaTildeStencils_[cell.id()] = gradGammaTilde_.grid.fluidCells().kNearestNeighbourSearch(cell.centroid(), 9);

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

    kappaMatrices_.resize(gamma_.size());
    for(const Cell &cell: gamma_.grid.fluidCells())
    {
        const size_t stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();

        A.resize(stencilSize, 5);

        int i = 0;
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar sSqr = (nb.cell().centroid() - cell.centroid()).magSqr();
            Scalar dx = nb.cell().centroid().x - cell.centroid().x;
            Scalar dy = nb.cell().centroid().y - cell.centroid().y;

            for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
            {
                if(ibObj.cells().isInGroup(nb.cell()))
                {
                    Point2D xc = ibObj.firstIntersect(cell.centroid(), nb.cell().centroid()).first;
                    sSqr = (xc - cell.centroid()).magSqr();
                    dx = xc.x - cell.centroid().x;
                    dy = xc.y - cell.centroid().y;
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
            Scalar sSqr = (dg.cell().centroid() - cell.centroid()).magSqr();
            Scalar dx = dg.cell().centroid().x - cell.centroid().x;
            Scalar dy = dg.cell().centroid().y - cell.centroid().y;

            for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
            {
                if(ibObj.cells().isInGroup(dg.cell()))
                {
                    Point2D xc = ibObj.firstIntersect(cell.centroid(), dg.cell().centroid()).first;
                    sSqr = (xc - cell.centroid()).magSqr();
                    dx = xc.x - cell.centroid().x;
                    dy = xc.y - cell.centroid().y;
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
            const Scalar sSqr = (bd.face().centroid() - cell.centroid()).magSqr();
            const Scalar dx = bd.face().centroid().x - cell.centroid().x;
            const Scalar dy = bd.face().centroid().y - cell.centroid().y;

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

//- Protected methods

void Celeste::computeGradGammaTilde()
{
    gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    // interpolateFaces(gammaTilde_);
    gradGammaTilde_.fill(Vector2D(0., 0.));

    Matrix b(8, 1);

    for(const Cell &cell: gradGammaTilde_.grid.fluidCells())
    {
        const auto &nb = gradGammaTildeStencils_[cell.id()];

        b.resize(8, 1);
        int i = 0;

        for(const Cell &kCell: nb)
        {
            if(&cell == &kCell)
                continue;

            const Scalar sSqr = (kCell.centroid() - cell.centroid()).magSqr();
            b(i, 0) = (gammaTilde_[kCell.id()] - gammaTilde_[cell.id()])/sSqr;

            ++i;
        }

        b = gradGammaTildeMatrices_[cell.id()]*b;
        gradGammaTilde_[cell.id()] = Vector2D(b(0, 0), b(1, 0));
    }

    //interpolateFaces(gradGammaTilde_);
}

void Celeste::computeInterfaceNormals()
{
    for(const Cell &cell: n_.grid.fluidCells())
        n_[cell.id()] = gradGammaTilde_[cell.id()] == Vector2D(0., 0.) ? Vector2D(0., 0.) : -gradGammaTilde_[cell.id()].unitVec();

    //for(const Face &face: n_.grid.faces())
    //    n_.faces()[face.id()] = gradGammaTilde_.faces()[face.id()] == Vector2D(0., 0.) ? Vector2D(0., 0.) : gradGammaTilde_.faces()[face.id()].unitVec();
}

void Celeste::computeCurvature()
{
    Matrix b(8, 1);
    kappa_.fill(0.);

    for(const Cell &cell: kappa_.grid.fluidCells())
    {
        const size_t stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();

        for(int compNo = 0; compNo < 2; ++compNo)
        {
            int i = 0;
            b.resize(stencilSize, 1);

            for(const InteriorLink &nb: cell.neighbours())
            {
                Scalar sSqr = (nb.cell().centroid() - cell.centroid()).magSqr();
                Scalar n = n_[nb.cell().id()](compNo) - n_[cell.id()](compNo);

                for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
                {
                    if(ibObj.cells().isInGroup(nb.cell()))
                    {
                        Point2D xc = ibObj.firstIntersect(cell.centroid(), nb.cell().centroid()).first;
                        sSqr = (xc - cell.centroid()).magSqr();
                        n = computeContactLineNormal(gradGammaTilde_[cell.id()], ibObj.centroid() - xc, u_[cell.id()])(compNo) - n_[cell.id()](compNo);
                        break;
                    }
                }

                b(i, 0) = n/sSqr;

                ++i;
            }

            for(const DiagonalCellLink &dg: cell.diagonals())
            {
                Scalar sSqr = (dg.cell().centroid() - cell.centroid()).magSqr();
                Scalar n = n_[dg.cell().id()](compNo) - n_[cell.id()](compNo);

                for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
                {
                    if(ibObj.cells().isInGroup(dg.cell()))
                    {
                        Point2D xc = ibObj.firstIntersect(cell.centroid(), dg.cell().centroid()).first;
                        sSqr = (xc - cell.centroid()).magSqr();
                        n = computeContactLineNormal(gradGammaTilde_[cell.id()], ibObj.centroid() - xc, u_[cell.id()])(compNo) - n_[cell.id()](compNo);
                        break;
                    }
                }

                b(i, 0) = n/sSqr;

                ++i;
            }

            for(const BoundaryLink &bd: cell.boundaries())
            {
                const Scalar sSqr = (bd.face().centroid() - cell.centroid()).magSqr();
                b(i, 0) = (n_.faces()[bd.face().id()](compNo) - n_[cell.id()](compNo))/sSqr;

                ++i;
            }

            b = kappaMatrices_[cell.id()]*b;

            kappa_[cell.id()] += b(compNo, 0);
        }
    }

    //weightCurvatures();
    interpolateFaces(kappa_);
}

void Celeste::weightCurvatures()
{
    auto pow8 = [](Scalar x){ return x*x*x*x*x*x*x*x; };

    for(const Cell &cell: wGamma_.grid.fluidCells())
        wGamma_[cell.id()] = pow8(1. - 2*fabs(0.5 - gammaTilde_[cell.id()]));

    kappa_.savePreviousIteration();

    for(const Cell &cell: kappa_.grid.fluidCells())
    {
        Scalar sumW = wGamma_[cell.id()], sumKappaW = wGamma_[cell.id()]*kappa_.prevIter()[cell.id()];

        for(const InteriorLink &nb: cell.neighbours())
        {
            sumW += wGamma_[nb.cell().id()];
            sumKappaW += kappa_.prevIter()[nb.cell().id()]*wGamma_[nb.cell().id()];
        }

        for(const DiagonalCellLink &dg: cell.diagonals())
        {
            sumW += wGamma_[dg.cell().id()];
            sumKappaW += kappa_.prevIter()[dg.cell().id()]*wGamma_[dg.cell().id()];
        }

        kappa_[cell.id()] = fabs(sumW) < 1e-15 ? 0. : sumKappaW/sumW;
    }

//    kappa_.savePreviousIteration();

//    for(const Cell &cell: kappa_.grid.fluidCells())
//    {
//        Scalar sumW = wGamma_[cell.id()], sumKappaW = wGamma_[cell.id()]*kappa_.prevIter()[cell.id()];

//        for(const InteriorLink &nb: cell.neighbours())
//        {
//            const Scalar mQ = pow8(dot(n_[cell.id()], nb.rCellVec().unitVec()));

//            sumW += wGamma_[nb.cell().id()]*mQ;
//            sumKappaW += kappa_.prevIter()[nb.cell().id()]*wGamma_[nb.cell().id()]*mQ;
//        }

//        for(const DiagonalCellLink &dg: cell.diagonals())
//        {
//            const Scalar mQ = pow8(dot(n_[cell.id()], dg.rCellVec().unitVec()));

//            sumW += wGamma_[dg.cell().id()]*mQ;
//            sumKappaW += kappa_.prevIter()[dg.cell().id()]*wGamma_[dg.cell().id()]*mQ;
//        }

//        Scalar kappa = sumKappaW/sumW;
//        kappa_[cell.id()] = isnan(kappa) ? 0. : kappa;
//    }
}
