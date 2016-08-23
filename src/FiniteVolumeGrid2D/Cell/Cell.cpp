#include "Cell.h"
#include "Polygon.h"
#include "Exception.h"

Cell::Cell(const std::vector<Label> &nodeIds, const std::vector<Node> &nodes)
{   
    isActive_ = true;
    isFluidCell_ = true;

    for (Label id: nodeIds)
        nodes_.push_back(std::cref(nodes[id]));

    std::vector<Point2D> vertices;

    for(const Node& node: nodes_)
        vertices.push_back(node);

    cellShape_ = Polygon(vertices);

    volume_ = cellShape_.area();
    centroid_ = cellShape_.centroid();

    if(volume_ < 0.)
        throw Exception("Cell", "Cell", "faces are not oriented in a counter-clockwise manner.");
}

void Cell::addDiagonalLink(const Cell &cell)
{
    diagonalLinks_.push_back(DiagonalCellLink(*this, cell));
}

bool Cell::isInCell(const Point2D &point) const
{
    return cellShape_.isInside(point);
}

//- Private methods

void Cell::addBoundaryLink(const Face& face)
{
    if(!face.isBoundary())
        throw Exception("Cell", "addBoundaryLink", "cannot add a boundary link to a non-boundary face.");

    boundaryLinks_.push_back(BoundaryLink(*this, face));
}

void Cell::addInteriorLink(const Face& face, const Cell& cell)
{
    if(!face.isInterior())
        throw Exception("Cell", "addInteriorLink", "cannot add an interior link to a non-interior face.");

    interiorLinks_.push_back(InteriorLink(*this, face, cell));
}

//- External functions
bool cellsShareFace(const Cell& cellA, const Cell& cellB)
{
    for(const InteriorLink& nb: cellA.neighbours())
    {
        if(nb.cell().id() == cellB.id() || &nb.cell() == &cellB)
            return true;
    }

    return false;
}
