#include "ContactLineGhostCellStencil.h"
#include "GhostCellImmersedBoundaryObject.h"

ContactLineGhostCellStencil::ContactLineGhostCellStencil(const Cell &cell,
                                                         const GhostCellImmersedBoundaryObject &ibObj,
                                                         const ScalarFiniteVolumeField &gamma,
                                                         Scalar theta)
        :
        GhostCellStencil(cell, ibObj, gamma.grid())
{
    if(bpGrad(gamma).magSqr() < 1e-8)
        return;

    Vector2D m = -bpGrad(gamma).unitVec();
    Vector2D nw = ibObj.nearestEdgeNormal(bp_);
    Vector2D cl = cross(m, nw) < 0. ? nw.rotate(M_PI_2 - theta) : nw.rotate(theta - M_PI_2);
    Line2D d(cell.centroid(), cl.normalVec());

    bool validStencil = false;

    int nCandidateNodes = 4;
    while(!validStencil)
        for(const Node &node: gamma.grid().findNearestNodes(bp_, nCandidateNodes++))
        {
            auto cells = node.cells();

            if (ibObj.noneInIb(cells.begin(), cells.end()))
            {
//                for(const Cell& ncell: cells)
//                {
//                    if((cell.centroid() - Vector2D(0.39, 0.99)).mag() < 1e-8)
//                        std::cout << ncell.centroid() << " \n";
//                }

                Point2D xc = d.nearestPoint(node);
                if (node.cellBounded(xc))
                {
                    cells_ = cells;
                    cells_.push_back(cell_);
                    ip_ = xc;
                    validStencil = true;
                    break;
                }
            }
        }

    Point2D x[] = {
            cells_[0].get().centroid(),
            cells_[1].get().centroid(),
            cells_[2].get().centroid(),
            cells_[3].get().centroid()
    };

    A_ = Matrix(4, 4, {
            x[0].x * x[0].y, x[0].x, x[0].y, 1.,
            x[1].x * x[1].y, x[1].x, x[1].y, 1.,
            x[2].x * x[2].y, x[2].x, x[2].y, 1.,
            x[3].x * x[3].y, x[3].x, x[3].y, 1.
    }).invert();

    neumannCoeffs_ = Matrix(1, 4, {ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * A_;
    neumannCoeffs_.push_back(-1.);
//
//    for(Scalar c: neumannCoeffs())
//        std::cout << c << " ";
//    std::cout << std::endl;
//
//        for (const GhostCellStencil &st: gcIbObj->stencils())
//        {
//            gammaTilde(st.cell()) = 0.;
//
//            //- Produce characteristic lines d1 and d2
//            Vector2D nb = -gcIbObj->nearestEdgeNormal(st.boundaryPoint());
//            Line2D d1(st.cell().centroid(), nb.rotate(M_PI_2 - theta).normalVec());
//            Line2D d2(st.cell().centroid(), nb.rotate(theta - M_PI_2).normalVec());
//
//            //- Track whether the line has been constructed
//            bool isConstructed[2] = {false, false};
//            std::vector<Ref<const Cell>> dCells[2];
//            Point2D ip[2];
//
//            const CellZone &fluid = grid_->cellZone("fluid");
//
//            //- Look for stencil candidates
//            int nNodes = 1;
//            while (!isConstructed[0] || !isConstructed[1])
//            {
//                //- Find candidate nodes
//                for (const Node &node: grid_->findNearestNodes(st.boundaryPoint(), nNodes++))
//                {
//                    auto cells = node.cells();
//
//                    if (fluid.isInGroup(cells.begin(), cells.end()))
//                    {
//                        Polygon bilinearStencil = Polygon({
//                                                                  cells[0].get().centroid(),
//                                                                  cells[1].get().centroid(),
//                                                                  cells[2].get().centroid(),
//                                                                  cells[3].get().centroid()
//                                                          }).convexHull();
//
//                        if (!isConstructed[0])
//                        {
//                            ip[0] = d1.nearestPoint(bilinearStencil.centroid());
//                            if (bilinearStencil.isInside(ip[0]))
//                            {
//                                dCells[0] = cells;
//                                isConstructed[0] = true;
//                            }
//                        }
//
//                        if (!isConstructed[1])
//                        {
//                            ip[1] = d2.nearestPoint(bilinearStencil.centroid());
//                            if (bilinearStencil.isInside(ip[1]))
//                            {
//                                dCells[1] = cells;
//                                isConstructed[1] = true;
//                            }
//                        }
//                    }
//                }
//            }
//
//            //- Stencils now constructed
//            for (int i = 0; i < 2; ++i)
//            {
//                Point2D x[] = {
//                        dCells[i][0].get().centroid(),
//                        dCells[i][1].get().centroid(),
//                        dCells[i][2].get().centroid(),
//                        dCells[i][3].get().centroid()
//                };
//
//                Matrix A(4, 4, {
//                        x[0].x * x[0].y, x[0].x, x[0].y, 1.,
//                        x[1].x * x[1].y, x[1].x, x[1].y, 1.,
//                        x[2].x * x[2].y, x[2].x, x[2].y, 1.,
//                        x[3].x * x[3].y, x[3].x, x[3].y, 1.
//                });
//
//                Matrix b(4, 1, {
//                        gammaTilde(dCells[i][0]),
//                        gammaTilde(dCells[i][1]),
//                        gammaTilde(dCells[i][2]),
//                        gammaTilde(dCells[i][3])
//                });
//
//                gammaTilde(st.cell()) += (Matrix(1, 4, {ip[i].x * ip[i].y, ip[i].x, ip[i].y, 1.}) * A.solve(b))(0, 0);
//                gammaTilde(st.cell()) = clamp(gammaTilde(st.cell()), 0., 1.);
//            }
//        }
}