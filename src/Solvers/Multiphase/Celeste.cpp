#include "Celeste.h"

Celeste::Celeste(const Input &input,
                 const ScalarFiniteVolumeField &gamma,
                 const VectorFiniteVolumeField &u,
                 std::map<std::string, ScalarFiniteVolumeField> &scalarFields,
                 std::map<std::string, VectorFiniteVolumeField> &vectorFields)
    :
      ContinuumSurfaceForce(input, gamma, u, scalarFields, vectorFields),
      wGamma_(gamma.grid, "wGamma")
{
    constructMatrices();
}

VectorFiniteVolumeField Celeste::compute()
{
    return ContinuumSurfaceForce::compute();
}

void Celeste::constructMatrices()
{
    Matrix A(8, 5);
    matrices_.resize(gamma_.size());

    for(const Cell &cell: gamma_.grid.fluidCells())
    {
        const size_t stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();

        A.resize(stencilSize, 5);

        int i = 0;
        for(const InteriorLink &nb: cell.neighbours())
        {
            const Scalar sSqr = (nb.cell().centroid() - cell.centroid()).magSqr();
            const Scalar dx = nb.cell().centroid().x - cell.centroid().x;
            const Scalar dy = nb.cell().centroid().y - cell.centroid().y;

            A(i, 0) = dx/sSqr;
            A(i, 1) = dy/sSqr;
            A(i, 2) = dx*dx/(2.*sSqr);
            A(i, 3) = dy*dy/(2.*sSqr);
            A(i, 4) = dx*dy/sSqr;

            ++i;
        }

        for(const DiagonalCellLink &dg: cell.diagonals())
        {
            const Scalar sSqr = (dg.cell().centroid() - cell.centroid()).magSqr();
            const Scalar dx = dg.cell().centroid().x - cell.centroid().x;
            const Scalar dy = dg.cell().centroid().y - cell.centroid().y;

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

        matrices_[cell.id()] = inverse(transpose(A)*A)*transpose(A);
    }
}

//- Protected methods

void Celeste::computeGradGammaTilde()
{
    gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    gradGammaTilde_.fill(Vector2D(0., 0.));

    Matrix b(8, 1);

    for(const Cell &cell: gradGammaTilde_.grid.fluidCells())
    {
        const size_t stencilSize = cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size();

            int i = 0;
            b.resize(stencilSize, 1);

            for(const InteriorLink &nb: cell.neighbours())
            {
                const Scalar sSqr = (nb.cell().centroid() - cell.centroid()).magSqr();
                b(i, 0) = (gammaTilde_[nb.cell().id()] - gammaTilde_[cell.id()])/sSqr;

                ++i;
            }

            for(const DiagonalCellLink &dg: cell.diagonals())
            {
                const Scalar sSqr = (dg.cell().centroid() - cell.centroid()).magSqr();
                b(i, 0) = (gammaTilde_[dg.cell().id()] - gammaTilde_[cell.id()])/sSqr;

                ++i;
            }

            for(const BoundaryLink &bd: cell.boundaries())
            {
                const Scalar sSqr = (bd.face().centroid() - cell.centroid()).magSqr();
                b(i, 0) = (gammaTilde_.faces()[bd.face().id()] - gammaTilde_[cell.id()])/sSqr;

                ++i;
            }

            b = matrices_[cell.id()]*b;
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
                const Scalar sSqr = (nb.cell().centroid() - cell.centroid()).magSqr();
                b(i, 0) = (n_[nb.cell().id()](compNo) - n_[cell.id()](compNo))/sSqr;

                ++i;
            }

            for(const DiagonalCellLink &dg: cell.diagonals())
            {
                const Scalar sSqr = (dg.cell().centroid() - cell.centroid()).magSqr();
                b(i, 0) = (n_[dg.cell().id()](compNo) - n_[cell.id()](compNo))/sSqr;

                ++i;
            }

            for(const BoundaryLink &bd: cell.boundaries())
            {
                const Scalar sSqr = (bd.face().centroid() - cell.centroid()).magSqr();
                b(i, 0) = (n_.faces()[bd.face().id()](compNo) - n_[cell.id()](compNo))/sSqr;

                ++i;
            }

            b = matrices_[cell.id()]*b;

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
