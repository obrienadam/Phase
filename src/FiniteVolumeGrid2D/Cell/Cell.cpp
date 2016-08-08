#include "Cell.h"
#include "Polygon.h"
#include "Exception.h"

Cell::Cell(const std::vector<Label> &nodeIds, const std::vector<Node> &nodes)
{   
    isActive_ = true;

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
