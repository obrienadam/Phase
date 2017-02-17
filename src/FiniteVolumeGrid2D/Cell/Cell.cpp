#include "Cell.h"
#include "Polygon.h"
#include "Exception.h"
#include "FiniteVolumeGrid2D.h"

Cell::Cell(const std::vector<Label> &nodeIds, const FiniteVolumeGrid2D &grid)
    :
      nodes_(grid.nodes()),
      nodeIds_(nodeIds)
{
    std::vector<Point2D> vertices;

    for(Label id: nodeIds_)
        vertices.push_back(nodes_[id]);

    cellShape_ = Polygon(vertices);

    volume_ = cellShape_.area();
    centroid_ = cellShape_.centroid();

    if(volume_ < 0.)
        throw Exception("Cell", "Cell", "faces are not oriented in a counter-clockwise manner.");

    id_ = grid.cells().size();
}

void Cell::addDiagonalLink(const Cell &cell)
{
    diagonalLinks_.push_back(DiagonalCellLink(*this, cell));
}

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

const std::vector< Ref<const Node> > Cell::nodes() const
{
    using namespace std;

    vector< Ref<const Node> > nodes;
    for(Label id: nodeIds_)
        nodes.push_back(cref(nodes_[id]));

    return nodes;
}

bool Cell::isInCell(const Point2D &point) const
{
    return cellShape_.isInside(point);
}

//- Private methods

//- External functions
bool cellsShareFace(const Cell& cellA, const Cell& cellB)
{
    for(const InteriorLink& nb: cellA.neighbours())
    {
        if(nb.cell().id() == cellB.id())
            return true;
    }

    return false;
}
