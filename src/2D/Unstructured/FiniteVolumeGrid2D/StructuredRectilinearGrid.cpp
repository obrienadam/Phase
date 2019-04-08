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

    std::vector<Point2D> xDimRefinements, yDimRefinements;

    //- Get children refinements, this is separated for compatibility reasons
    printf("Applying refinements...\n");
    auto xRefinementInput = input.caseInput().get_child_optional("Grid.refineX");

    if(xRefinementInput)
        for(const auto& entry: xRefinementInput.get())
            xDimRefinements.emplace_back(entry.first);

    auto yRefinementInput = input.caseInput().get_child_optional("Grid.refineY");

    if(yRefinementInput)
        for(const auto& entry: yRefinementInput.get())
            yDimRefinements.emplace_back(entry.first);

    Point2D origin = input.caseInput().get<std::string>("Grid.origin", "(0,0)");

    init(width, height, nCellsX, nCellsY, convertToMeters, xDimRefinements, yDimRefinements, origin);
}

StructuredRectilinearGrid::StructuredRectilinearGrid(const std::vector<Scalar> &xcoords, const std::vector<Scalar> &ycoords)
{
    printf("Initiating grid from coordinates...\n");
    init(xcoords, ycoords);
}

void StructuredRectilinearGrid::init(Scalar width, Scalar height,
                                     Size nCellsX, Size nCellsY,
                                     Scalar convertToMeters,
                                     const std::vector<Point2D> &xDimRefinements,
                                     const std::vector<Point2D> &yDimRefinements,
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
        xDims.emplace_back(i * hx0);

    for (Label j = 0; j < nCellsY + 1; ++j)
        yDims.emplace_back(j * hy0);

    for (const auto &xDimRefinement: xDimRefinements)
        xDims = refineDims(xDimRefinement.x, xDimRefinement.y, xDims);

    for (const auto &yDimRefinement: yDimRefinements)
        yDims = refineDims(yDimRefinement.x, yDimRefinement.y, yDims);

    Size nNodesX = xDims.size();
    Size nNodesY = yDims.size();
    nCellsX_ = nNodesX - 1;
    nCellsY_ = nNodesY - 1;

    std::vector<Point2D> nodes;
    nodes.reserve(nNodesX * nNodesY);
    for (Label j = 0; j < nNodesY; ++j)
        for (Label i = 0; i < nNodesX; ++i)
            nodes.emplace_back(xDims[i], yDims[j]);

    //- Create cells
    std::vector<Label> elemInds(1, 0), elems;
    elemInds.reserve(nCellsX_ * nCellsY_);
    elems.reserve(4 * nCellsX_ * nCellsY_);

    for (Label j = 0; j < nCellsY_; ++j)
        for (Label i = 0; i < nCellsX_; ++i)
        {
            elemInds.emplace_back(elemInds.back() + 4);
            elems.emplace_back(j * nNodesX + i);
            elems.emplace_back(j * nNodesX + i + 1);
            elems.emplace_back((j + 1) * nNodesX + i + 1);
            elems.emplace_back((j + 1) * nNodesX + i);
        }

    FiniteVolumeGrid2D::init(nodes, elemInds, elems, origin);
    initPatches();
}

void StructuredRectilinearGrid::init(std::vector<Scalar> xcoords,
                                     std::vector<Scalar> ycoords)
{
    if(xcoords.size() < 2 || ycoords.size() < 2)
        throw Exception("StructuredRectilinearGrid", "init", "invalid grid coordinates.");

    std::sort(xcoords.begin(), xcoords.end());
    std::sort(ycoords.begin(), ycoords.end());

    width_ = xcoords.back() - xcoords.front();
    height_ = ycoords.back() - ycoords.front();
    nCellsX_ = xcoords.size() - 1;
    nCellsY_ = ycoords.size() - 1;

    std::vector<Point2D> nodes;
    nodes.reserve(xcoords.size() * ycoords.size());

    for(size_t j = 0; j < ycoords.size(); ++j)
        for(size_t i = 0; i < xcoords.size(); ++i)
            nodes.emplace_back(xcoords[i], ycoords[j]);

    //- Create cells
    std::vector<Label> elemInds(1, 0), elems;
    elemInds.reserve(nCellsX_ * nCellsY_);
    elems.reserve(4 * nCellsX_ * nCellsY_);

    for (Label j = 0; j < nCellsY_; ++j)
        for (Label i = 0; i < nCellsX_; ++i)
        {
            elemInds.emplace_back(elemInds.back() + 4);
            elems.emplace_back(j * xcoords.size() + i);
            elems.emplace_back(j * xcoords.size() + i + 1);
            elems.emplace_back((j + 1) * xcoords.size() + i + 1);
            elems.emplace_back((j + 1) * xcoords.size() + i);
        }

    FiniteVolumeGrid2D::init(nodes, elemInds, elems, Point2D(0., 0.));
    initPatches();
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

std::vector<Scalar> StructuredRectilinearGrid::refineDims(Scalar start, Scalar end, const std::vector<Scalar> &dims)
{
    std::vector<Scalar> newDims;

    newDims.reserve(dims.size());
    for (int i = 0; i < dims.size() - 1; ++i)
    {
        Scalar x = dims[i];
        newDims.emplace_back(x);

        if (x >= start && x < end)
            newDims.emplace_back(x + (dims[i + 1] - x) / 2.);
    }

    newDims.emplace_back(dims.back());

    return newDims;
}

void StructuredRectilinearGrid::initPatches()
{
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
