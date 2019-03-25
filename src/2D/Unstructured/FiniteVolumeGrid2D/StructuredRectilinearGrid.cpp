#include "StructuredRectilinearGrid.h"

StructuredRectilinearGrid::StructuredRectilinearGrid(const Input &input)
        :
        FiniteVolumeGrid2D()
{
    Scalar width = input.caseInput().get<Scalar>("Grid.width");
    Scalar height = input.caseInput().get<Scalar>("Grid.height");
    size_t nCellsX = input.caseInput().get<size_t>("Grid.nCellsX");
    size_t nCellsY = input.caseInput().get<size_t>("Grid.nCellsY");
    Scalar convertToMeters = input.caseInput().get<Scalar>("Grid.convertToMeters", 1.);

    std::vector<std::pair<Scalar, Scalar>> xDimRefinements, yDimRefinements;

    Vector2D tmp = input.caseInput().get<std::string>("Grid.refineX", "(0,0)");
    xDimRefinements.push_back(std::make_pair(tmp.x, tmp.y));

    tmp = input.caseInput().get<std::string>("Grid.refineY", "(0,0)");
    yDimRefinements.push_back(std::make_pair(tmp.x, tmp.y));

    Point2D origin = input.caseInput().get<std::string>("Grid.origin", "(0,0)");

    init(width, height, nCellsX, nCellsY, convertToMeters, xDimRefinements, yDimRefinements, origin);
}

void StructuredRectilinearGrid::init(Scalar width, Scalar height,
                                     Size nCellsX, Size nCellsY,
                                     Scalar convertToMeters,
                                     const std::vector<std::pair<Scalar, Scalar> > &xDimRefinements,
                                     const std::vector<std::pair<Scalar, Scalar> > &yDimRefinements,
                                     const Point2D &origin)
{
    width *= convertToMeters;
    height *= convertToMeters;

    width_ = width;
    height_ = height;

    Scalar hx0 = width_ / nCellsX;
    Scalar hy0 = height_ / nCellsY;

    //- Create nodes
    std::vector<Scalar> xDims, yDims;

    for (Label i = 0; i < nCellsX + 1; ++i)
        xDims.push_back(i * hx0);

    for (Label j = 0; j < nCellsY + 1; ++j)
        yDims.push_back(j * hy0);

    for (const auto &xDimRefinement: xDimRefinements)
        refineDims(xDimRefinement.first, xDimRefinement.second, xDims);

    for (const auto &yDimRefinement: yDimRefinements)
        refineDims(yDimRefinement.first, yDimRefinement.second, yDims);

    Size nNodesX = xDims.size();
    Size nNodesY = yDims.size();
    nCellsX_ = nNodesX - 1;
    nCellsY_ = nNodesY - 1;

    std::vector<Point2D> nodes;
    for (Label j = 0; j < nNodesY; ++j)
        for (Label i = 0; i < nNodesX; ++i)
            nodes.push_back(Point2D(xDims[i], yDims[j]));

    //- Create cells
    std::vector<Label> elemInds(1, 0), elems;

    for (Label j = 0; j < nCellsY_; ++j)
        for (Label i = 0; i < nCellsX_; ++i)
        {
            elemInds.push_back(elemInds.back() + 4);
            elems.push_back(j * nNodesX + i);
            elems.push_back(j * nNodesX + i + 1);
            elems.push_back((j + 1) * nNodesX + i + 1);
            elems.push_back((j + 1) * nNodesX + i);
        }

    FiniteVolumeGrid2D::init(nodes, elemInds, elems, origin);

    //- Construct default patches
    std::vector<Label> xm, xp, ym, yp;

    //- x patches
    for (Label j = 0; j < nCellsY_; ++j)
    {
        xm.push_back(
                findFace(node(0, j).id(), node(0, j + 1).id())
        );
        xp.push_back(
                findFace(node(nCellsX_, j).id(), node(nCellsX_, j + 1).id())
        );
    }
    createPatch("x-", xm);
    createPatch("x+", xp);

    //- y patches
    for (Label i = 0; i < nCellsX_; ++i)
    {
        ym.push_back(
                findFace(node(i, 0).id(), node(i + 1, 0).id())
        );
        yp.push_back(
                findFace(node(i, nCellsY_).id(), node(i + 1, nCellsY_).id())
        );
    }
    createPatch("y-", ym);
    createPatch("y+", yp);
}


Cell &StructuredRectilinearGrid::operator()(Label i, Label j)
{
    return cells_[nCellsX_ * j + i];
}

const Cell &StructuredRectilinearGrid::operator()(Label i, Label j) const
{
    return cells_[nCellsX_ * j + i];
}

const Node &StructuredRectilinearGrid::node(Label i, Label j) const
{
    return nodes_[(nCellsX_ + 1) * j + i];
}

Scalar StructuredRectilinearGrid::h() const
{
    if (!isEquidistant())
        throw Exception("StructuredRectilinearGrid", "h", "not an equidistant grid.");

    return width_ / nCellsX_;
}

void StructuredRectilinearGrid::refineDims(Scalar start, Scalar end, std::vector<Scalar> &dims)
{
    std::vector<Scalar> newDims;

    for (int i = 0; i < dims.size() - 1; ++i)
    {
        Scalar x = dims[i];

        newDims.push_back(x);

        if (x >= start && x < end)
            newDims.push_back(x + (dims[i + 1] - x) / 2.);
    }

    newDims.push_back(dims.back());

    dims = std::move(newDims);
}
