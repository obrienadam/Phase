#include "Celeste.h"
#include "Matrix.h"

Celeste::Celeste(const Input &input,
                 const ScalarFiniteVolumeField &gamma,
                 const VectorFiniteVolumeField &u,
                 std::map<std::string, ScalarFiniteVolumeField> &scalarFields,
                 std::map<std::string, VectorFiniteVolumeField> &vectorFields)
    :
      ContinuumSurfaceForce(input, gamma, u, scalarFields, vectorFields),
      wGamma_(gamma.grid, "wGamma")
{

}

VectorFiniteVolumeField Celeste::compute()
{
    return ContinuumSurfaceForce::compute();
}

//- Protected methods

void Celeste::computeGradGammaTilde()
{
    gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    gradGammaTilde_.fill(Vector2D(0., 0.));

    Matrix A(8, 5), b(8, 1);

    for(const Cell &cell: gammaTilde_.grid.fluidCells())
    {
        auto nb = gammaTilde_.grid.activeCells().kNearestNeighbourSearch(cell.centroid(), 9);

        int i = 0;
        for(const Cell &kCell: nb)
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

            b(i, 0) = (gammaTilde_[kCell.id()] - gammaTilde_[cell.id()])/sSqr;

            ++i;
        }

        A.solve(b);

        gradGammaTilde_[cell.id()] = Vector2D(b(0, 0), b(1, 0));
    }

    interpolateFaces(gradGammaTilde_);
}

void Celeste::computeInterfaceNormals()
{
    ContinuumSurfaceForce::computeInterfaceNormals();
}

void Celeste::computeCurvature()
{
    kappa_.fill(0.);

    Matrix A(8, 5), b(8, 1);

    for(const Cell &cell: kappa_.grid.fluidCells())
    {
        auto nb = kappa_.grid.activeCells().kNearestNeighbourSearch(cell.centroid(), 9);

        A.resize(8 + cell.boundaries().size(), 5);
        b.resize(8 + cell.boundaries().size(), 1);

        for(int compNo = 0; compNo < 2; ++compNo)
        {
            int i = 0;
            for(const Cell &kCell: nb)
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

                b(i, 0) = (n_[kCell.id()](compNo) - n_[cell.id()](compNo))/sSqr;

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

                b(i, 0) = (n_.faces()[bd.face().id()](compNo) - n_[cell.id()](compNo))/sSqr;

                ++i;
            }

            A.solve(b);

            kappa_[cell.id()] += b(compNo, 0);
        }
    }

    //weightCurvatures();
    interpolateFaces(kappa_);
}

void Celeste::weightCurvatures()
{
    auto pow8 = [](Scalar x){ return x*x*x*x*x*x*x*x; };

    for(const Cell &cell: wGamma_.grid.activeCells())
        wGamma_[cell.id()] = pow8(1. - 2*fabs(0.5 - gamma_[cell.id()]));

    kappa_.savePreviousIteration();

    for(const Cell &cell: kappa_.grid.fluidCells())
    {
        Scalar sumW = wGamma_[cell.id()], sumKappaW = wGamma_[cell.id()]*kappa_.prevIter()[cell.id()];

        for(const InteriorLink &nb: cell.neighbours())
        {
            sumW += wGamma_[nb.cell().id()];
            sumKappaW += kappa_.prevIter()[nb.cell().id()]*wGamma_[nb.cell().id()];
        }

        kappa_[cell.id()] = isnan(sumKappaW/sumW) ? 0. : sumKappaW/sumW;
    }

    kappa_.savePreviousIteration();

    for(const Cell &cell: kappa_.grid.fluidCells())
    {
        Scalar sumW = wGamma_[cell.id()], sumKappaW = wGamma_[cell.id()]*kappa_.prevIter()[cell.id()];

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Scalar mQ = pow8(dot(n_[cell.id()], nb.rCellVec().unitVec()));

            sumW += wGamma_[nb.cell().id()]*mQ;
            sumKappaW += kappa_.prevIter()[nb.cell().id()]*wGamma_[nb.cell().id()]*mQ;
        }

        kappa_[cell.id()] = isnan(sumKappaW/sumW) ? 0. : sumKappaW/sumW;
    }
}
