#include "StructuredRectilinearGrid.h"
#include "Exception.h"

StructuredRectilinearGrid::StructuredRectilinearGrid(Scalar width, Scalar height,
                                                     Size nCellsX, Size nCellsY,
                                                     Scalar convertToMeters,
                                                     const std::vector<std::pair<Scalar, Scalar> > &xDimRefinements,
                                                     const std::vector<std::pair<Scalar, Scalar> > &yDimRefinements)
    :
      FiniteVolumeGrid2D((nCellsX + 1)*(nCellsY + 1),
                         nCellsX*nCellsY,
                         (nCellsX + 1)*nCellsY + nCellsX*(nCellsY + 1))
{
    width *= convertToMeters;
    height *= convertToMeters;

    width_ = width;
    height_ = height;

    Scalar hx0 = width_/nCellsX;
    Scalar hy0 = height_/nCellsY;

    //- Create nodes
    std::vector<Scalar> xDims, yDims;

    for(Label i = 0; i < nCellsX + 1; ++i)
        xDims.push_back(i*hx0);

    for(Label j = 0; j < nCellsY + 1; ++j)
        yDims.push_back(j*hy0);

    for(const auto& xDimRefinement: xDimRefinements)
        refineDims(xDimRefinement.first, xDimRefinement.second, xDims);

    for(const auto& yDimRefinement: yDimRefinements)
        refineDims(yDimRefinement.first, yDimRefinement.second, yDims);

    Size nNodesX = xDims.size();
    Size nNodesY = yDims.size();
    nCellsX_ = nNodesX - 1;
    nCellsY_ = nNodesY - 1;

    std::vector<Point2D> nodes;
    for(Label j = 0; j < nNodesY; ++j)
        for(Label i = 0; i < nNodesX; ++i)
            nodes.push_back(Point2D(xDims[i], yDims[j]));

    //- Create cells
    std::vector<Label> elemInds(1, 0), elems;

    for(Label j = 0; j < nCellsY_; ++j)
        for(Label i = 0; i < nCellsX_; ++i)
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

void StructuredRectilinearGrid::refineDims(Scalar start, Scalar end, std::vector<Scalar> &dims)
{
    std::vector<Scalar> newDims;

    for(int i = 0; i < dims.size() - 1; ++i)
    {
        Scalar x = dims[i];

        newDims.push_back(x);

        if(x >= start && x < end)
            newDims.push_back(x + (dims[i + 1] - x)/2.);
    }

    newDims.push_back(dims.back());

    dims = std::move(newDims);
}
