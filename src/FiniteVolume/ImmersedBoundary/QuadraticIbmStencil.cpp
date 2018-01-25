#include "QuadraticIbmStencil.h"
#include "StaticMatrix.h"

QuadraticIbmStencil::QuadraticIbmStencil(const InteriorLink &link,
                                         const ImmersedBoundary &ib)
{
    auto ibObj = ib.ibObj(link.cell().centroid());

    if (ibObj)
    {
        const Cell &stCell = ib.grid()->globalActiveCells().nearestItem(
                2. * link.self().centroid() - link.cell().centroid()
        );

        auto ibObjStCell = ib.ibObj(stCell.centroid());

        if (ibObjStCell)
            initQuadraticCoeffs(stCell, link.self(), link.cell(), *ibObjStCell, *ibObj);
        else
            initQuadraticCoeffs(stCell, link.self(), link.cell(), *ibObj);
    }
    else
    {
        cells_ = {std::cref(link.cell())};
        coeffs_ = {1.};
    }
}

void QuadraticIbmStencil::initQuadraticCoeffs(const Cell &stCell,
                                              const Cell &cell,
                                              const Cell &ibCell,
                                              const ImmersedBoundaryObject &ibObj)
{
    Vector2D eta = (ibCell.centroid() - stCell.centroid()).unitVec();
    LineSegment2D ln = ibObj.intersectionLine(cell.centroid(), ibCell.centroid());

    Scalar e[] = {
            dot(stCell.centroid(), eta),
            dot(cell.centroid(), eta),
            dot(ln.ptB(), eta)
    };

    try
    {
        auto A = inverse(StaticMatrix<3, 3>({
                                                    e[0] * e[0], e[0], 1.,
                                                    e[1] * e[1], e[1], 1.,
                                                    e[2] * e[2], e[2], 1.
                                            }));

        Scalar eb = dot(ibCell.centroid(), eta);
        auto c = StaticMatrix<1, 3>({eb * eb, eb, 1.}) * A;

        cells_ = {
                std::cref(stCell),
                std::cref(cell)
        };

        coeffs_ = {
                c(0, 0),
                c(0, 1)
        };

        src_ = c(0, 2) * ibObj.velocity(ln.ptB());
    }
    catch (const Exception &e)
    {
        initLinearCoeffs(stCell, ibCell, ibObj);
    }
}

void QuadraticIbmStencil::initQuadraticCoeffs(const Cell &stCell,
                                              const Cell &cell,
                                              const Cell &ibCell,
                                              const ImmersedBoundaryObject &ibObjL,
                                              const ImmersedBoundaryObject &ibObjR)
{
    Vector2D eta = (ibCell.centroid() - stCell.centroid()).unitVec();
    LineSegment2D lnSt = ibObjL.intersectionLine(cell.centroid(), stCell.centroid());
    LineSegment2D lnIb = ibObjR.intersectionLine(cell.centroid(), ibCell.centroid());

    Scalar e[] = {
            dot(lnSt.ptB(), eta),
            dot(cell.centroid(), eta),
            dot(lnIb.ptB(), eta)
    };

    try
    {
        auto A = inverse(StaticMatrix<3, 3>({
                                                    e[0] * e[0], e[0], 1.,
                                                    e[1] * e[1], e[1], 1.,
                                                    e[2] * e[2], e[2], 1.
                                            }));

        Scalar eb = dot(ibCell.centroid(), eta);
        auto c = StaticMatrix<1, 3>({eb * eb, eb, 1.}) * A;

        cells_ = {
                std::cref(cell)
        };

        coeffs_ = {
                c(0, 1)
        };

        src_ = c(0, 0) * ibObjL.velocity(lnSt.ptB()) + c(0, 2) * ibObjR.velocity(lnIb.ptB());
    }
    catch (const Exception &e)
    {
        initLinearCoeffs(stCell, ibCell, ibObjL, ibObjR);
    }
}

void QuadraticIbmStencil::initLinearCoeffs(const Cell &stCell, const Cell &ibCell, const ImmersedBoundaryObject &ibObj)
{
    Vector2D eta = (ibCell.centroid() - stCell.centroid()).unitVec();
    LineSegment2D ln = ibObj.intersectionLine(stCell.centroid(), ibCell.centroid());

    Scalar e[] = {
            dot(stCell.centroid(), eta),
            dot(ln.ptB(), eta)
    };

    auto A = inverse(StaticMatrix<2, 2>({
                                                e[0], 1.,
                                                e[1], 1.
                                        }));

    Scalar eb = dot(ibCell.centroid(), eta);
    auto c = StaticMatrix<1, 2>({eb, 1.}) * A;

    cells_ = {
            stCell
    };

    coeffs_ = {
            c(0, 0)
    };

    src_ = c(0, 1) * ibObj.velocity(ln.ptB());
}

void QuadraticIbmStencil::initLinearCoeffs(const Cell &stCell,
                                           const Cell &ibCell,
                                           const ImmersedBoundaryObject &ibObjL,
                                           const ImmersedBoundaryObject &ibObjR)
{
    Vector2D eta = (ibCell.centroid() - stCell.centroid()).unitVec();
    LineSegment2D lnSt = ibObjL.intersectionLine(ibCell.centroid(), stCell.centroid());
    LineSegment2D lnIb = ibObjR.intersectionLine(stCell.centroid(), ibCell.centroid());

    Scalar e[] = {
            dot(lnSt.ptB(), eta),
            dot(lnIb.ptB(), eta)
    };

    auto A = inverse(StaticMatrix<2, 2>({
                                                e[0], 1.,
                                                e[1], 1.
                                        }));

    Scalar eb = dot(ibCell.centroid(), eta);
    auto c = StaticMatrix<1, 2>({eb, 1.}) * A;

    coeffs_ = {};
    src_ = c(0, 0) * ibObjL.velocity(lnSt.ptB()) + c(0, 1) * ibObjR.velocity(lnIb.ptB());
}