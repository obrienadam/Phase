#include "StructuredRectilinearGrid.h"
#include "Exception.h"

StructuredRectilinearGrid::StructuredRectilinearGrid(Scalar width, Scalar height, Size nCellsX, Size nCellsY, Scalar convertToMeters)
    :
      FiniteVolumeGrid2D((nCellsX + 1)*(nCellsY + 1),
                         nCellsX*nCellsY,
                         (nCellsX + 1)*nCellsY + nCellsX*(nCellsY + 1))
{
    width *= convertToMeters;
    height *= convertToMeters;

    nCellsX_ = nCellsX;
    nCellsY_ = nCellsY;
    width_ = width;
    height_ = height;

    Size nNodesX = nCellsX + 1;
    Size nNodesY = nCellsY + 1;
    Scalar hx = width/nCellsX_;
    Scalar hy = height/nCellsY_;

    //- Create nodes
    std::vector<Point2D> nodes;

    for(Label j = 0; j < nNodesY; ++j)
        for(Label i = 0; i < nNodesX; ++i)
            nodes.push_back(Point2D(i*hx, j*hy));

    //- Create cells
    std::vector<Label> elemInds(1, 0), elems;

    for(Label j = 0; j < nNodesY - 1; ++j)
        for(Label i = 0; i < nNodesX - 1; ++i)
        {
            elemInds.push_back(elemInds.back() + 4);
            elems.push_back(j*nNodesX + i);
            elems.push_back(j*nNodesX + i + 1);
            elems.push_back((j + 1)*nNodesX + i + 1);
            elems.push_back((j + 1)*nNodesX + i);
        }

    init(nodes, elemInds, elems);

    //- Construct default patches
    std::vector<Label> xm, xp, ym, yp;

    //- x patches
    for(Label j = 0; j < nCellsY_; ++j)
    {
        xm.push_back(
                    findFace(node(0, j).id(), node(0, j + 1).id())
                    );
        xp.push_back(
                    findFace(node(nCellsX_, j).id(), node(nCellsX_, j + 1).id())
                    );
    }
    applyPatch("x-", xm);
    applyPatch("x+", xp);

    //- y patches
    for(Label i = 0; i < nCellsX_; ++i)
    {
        ym.push_back(
                    findFace(node(i, 0).id(), node(i + 1, 0).id())
                    );
        yp.push_back(
                    findFace(node(i, nCellsY_).id(), node(i + 1, nCellsY_).id())
                    );
    }
    applyPatch("y-", ym);
    applyPatch("y+", yp);
}

Cell& StructuredRectilinearGrid::operator()(Label i, Label j)
{
    if(i < 0 || i >= nCellsX_
            || j < 0 || j >= nCellsY_)
        throw Exception("StructuredRectilinearGrid", "operator()", "index is out of range.");

    return cells_[nCellsX_*j + i];
}

const Cell& StructuredRectilinearGrid::operator()(Label i, Label j) const
{
    if(i < 0 || i >= nCellsX_
            || j < 0 || j >= nCellsY_)
        throw Exception("StructuredRectilinearGrid", "operator()", "index is out of range.");

    return cells_[nCellsX_*j + i];
}

const Node& StructuredRectilinearGrid::node(Label i, Label j) const
{
    if(i < 0 || i > nCellsX_
            || j < 0 || j > nCellsY_)
        throw Exception("StructuredRectilinearGrid", "node", "index is out of range.");

    return nodes_[(nCellsX_ + 1)*j + i];
}
