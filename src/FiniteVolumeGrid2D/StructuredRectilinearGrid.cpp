#include "StructuredRectilinearGrid.h"
#include "Exception.h"

StructuredRectilinearGrid::StructuredRectilinearGrid(Scalar width, Scalar height, Size nCellsX, Size nCellsY)
    :
      FiniteVolumeGrid2D((nCellsX + 1)*(nCellsY + 1),
                         nCellsX*nCellsY,
                         (nCellsX + 1)*nCellsY + nCellsX*(nCellsY + 1))
{
    nCellsX_ = nCellsX;
    nCellsY_ = nCellsY;
    width_ = width;
    height_ = height;

    Size nNodesX = nCellsX + 1;
    Size nNodesY = nCellsY + 1;
    Scalar hx = width/nCellsX_;
    Scalar hy = height/nCellsY_;

    //- Create nodes
    for(Label j = 0; j < nNodesY; ++j)
        for(Label i = 0; i < nNodesX; ++i)
            addNode(Point2D(i*hx, j*hy));

    //- Create cells
    for(Label j = 0; j < nCellsY_; ++j)
        for(Label i = 0; i < nCellsX_; ++i)
        {
            std::vector<Label> nids = {
                node(i, j).id(),
                node(i + 1, j).id(),
                node(i + 1, j + 1).id(),
                node(i, j + 1).id(),
            };

            createCell(nids);
        }

    //- Construct default patches
    std::vector< Ref<Face> > xm, xp, ym, yp;

    // x- patch
    for(Label j = 0; j < nCellsY_; ++j)
    {
        xm.push_back(
                    std::ref(faces_[findFace(node(0, j).id(), node(0, j + 1).id())])
                    );
    }
    applyPatch("x-", xm);

    // x+ patch
    for(Label j = 0; j < nCellsY_; ++j)
    {
        xp.push_back(
                    std::ref(faces_[findFace(node(nCellsX_, j).id(), node(nCellsX_, j + 1).id())])
                    );
    }
    applyPatch("x+", xp);

    // y- patch
    for(Label i = 0; i < nCellsX_; ++i)
    {
        ym.push_back(
                    std::ref(faces_[findFace(node(i, 0).id(), node(i + 1, 0).id())])
                    );
    }
    applyPatch("y-", ym);

    // y+ patch
    for(Label i = 0; i < nCellsX_; ++i)
    {
        yp.push_back(
                    std::ref(faces_[findFace(node(i, nCellsY_).id(), node(i + 1, nCellsY_).id())])
                    );
    }
    applyPatch("y+", yp);

    initConnectivity();
    computeBoundingBox();
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
