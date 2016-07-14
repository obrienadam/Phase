#include "ImmersedBoundaryCeleste.h"
#include "Matrix.h"
#include "ImmersedBoundaryObject.h"

//- Public methods

ImmersedBoundaryCeleste::ImmersedBoundaryCeleste(const Input &input,
                                                 const ScalarFiniteVolumeField &gamma,
                                                 const VectorFiniteVolumeField &u,
                                                 std::map<std::string, ScalarFiniteVolumeField> &scalarFields,
                                                 std::map<std::string, VectorFiniteVolumeField> &vectorFields)
    :
      Celeste(input, gamma, u, scalarFields, vectorFields)
{

}

VectorFiniteVolumeField ImmersedBoundaryCeleste::compute(const std::vector<ImmersedBoundaryObject> &ibObjs)
{
    Celeste::computeGradGammaTilde();
    Celeste::computeInterfaceNormals();
    computeCurvature(ibObjs);

    VectorFiniteVolumeField ft(gamma_.grid, "ft");
    VectorFiniteVolumeField gradGamma = grad(gamma_);

    for(const Cell &cell: gamma_.grid.fluidCells())
    {
        ft[cell.id()] = sigma_*kappa_[cell.id()]*gradGamma[cell.id()];
        if(isnan(ft[cell.id()].magSqr()))
            ft[cell.id()] = Vector2D(0, 0);
    }

    return ft;
}

//- Protected methods

void ImmersedBoundaryCeleste::computeCurvature(const std::vector<ImmersedBoundaryObject> &ibObjs)
{
    kappa_.fill(0.);

    Matrix A(8, 5), b(8, 1);

    for(const Cell &cell: kappa_.grid.fluidCells())
    {
        const size_t stencilSize = 8;

        A.resize(stencilSize, 5);
        b.resize(stencilSize, 1);

        for(int compNo = 0; compNo < 2; ++compNo)
        {
            int i = 0;
            for(const InteriorLink &nb: cell.neighbours())
            {
                Scalar sSqr = (nb.cell().centroid() - cell.centroid()).magSqr();
                Scalar dx = nb.cell().centroid().x - cell.centroid().x;
                Scalar dy = nb.cell().centroid().y - cell.centroid().y;
                Scalar dn = n_[nb.cell().id()](compNo) - n_[cell.id()](compNo);

                for(const ImmersedBoundaryObject &ibObj: ibObjs)
                {
                    if(ibObj.cells().isInGroup(nb.cell()))
                    {
                        Point2D xc = ibObj.firstIntersect(cell.centroid(), nb.cell().centroid()).first;
                        sSqr = (xc - cell.centroid()).magSqr();
                        dx = xc.x - cell.centroid().x;
                        dy = xc.y - cell.centroid().y;
                        dn = computeContactLineNormal(gammaTilde_[cell.id()], ibObj.centroid() - xc, u_[cell.id()])(compNo) - n_[cell.id()](compNo);
                    }
                }

                A(i, 0) = dx/sSqr;
                A(i, 1) = dy/sSqr;
                A(i, 2) = dx*dx/(2.*sSqr);
                A(i, 3) = dy*dy/(2.*sSqr);
                A(i, 4) = dx*dy/sSqr;

                b(i, 0) = dn/sSqr;

                ++i;
            }

            for(const DiagonalCellLink &dg: cell.diagonals())
            {
                Scalar sSqr = (dg.cell().centroid() - cell.centroid()).magSqr();
                Scalar dx = dg.cell().centroid().x - cell.centroid().x;
                Scalar dy = dg.cell().centroid().y - cell.centroid().y;
                Scalar dn = n_[dg.cell().id()](compNo) - n_[cell.id()](compNo);

                for(const ImmersedBoundaryObject &ibObj: ibObjs)
                {
                    if(ibObj.cells().isInGroup(dg.cell()) || !dg.cell().isActive())
                    {
                        Point2D xc = ibObj.firstIntersect(cell.centroid(), dg.cell().centroid()).first;
                        sSqr = (xc - cell.centroid()).magSqr();
                        dx = xc.x - cell.centroid().x;
                        dy = xc.y - cell.centroid().y;
                        dn = computeContactLineNormal(gammaTilde_[cell.id()], ibObj.centroid() - xc, u_[cell.id()])(compNo) - n_[cell.id()](compNo);
                    }
                }

                A(i, 0) = dx/sSqr;
                A(i, 1) = dy/sSqr;
                A(i, 2) = dx*dx/(2.*sSqr);
                A(i, 3) = dy*dy/(2.*sSqr);
                A(i, 4) = dx*dy/sSqr;

                b(i, 0) = dn/sSqr;

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
