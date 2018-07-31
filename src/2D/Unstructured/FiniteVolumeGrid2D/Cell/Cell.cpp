#include "Geometry/Polygon.h"
#include "System/Exception.h"
#include "FiniteVolumeGrid2D/FiniteVolumeGrid2D.h"

#include "Cell.h"

Cell::Cell(const std::vector<Label> &nodeIds, const FiniteVolumeGrid2D &grid)
        :
        grid_(grid),
        nodeIds_(nodeIds)
{
    std::vector<Point2D> vertices;
    vertices.reserve(nodeIds.size());

    for (Label id: nodeIds_)
        vertices.push_back(grid.nodes()[id]);

    cellShape_ = Polygon(vertices.begin(), vertices.end());

    volume_ = cellShape_.area();
    centroid_ = cellShape_.centroid();

    if (volume_ < 0.)
        throw Exception("Cell", "Cell", "faces are not oriented in a counter-clockwise manner.");

    localId_ = grid.cells().size();
    globalId_ = localId_;
}

Scalar Cell::polarVolume() const
{
    Scalar volume = 0.;

    for (const InteriorLink &nb: interiorLinks_)
        volume += dot(nb.face().centroid(), nb.face().polarOutwardNorm(centroid_));

    for (const BoundaryLink &bd: boundaryLinks_)
        volume += dot(bd.face().centroid(), bd.face().polarOutwardNorm(centroid_));

//    Scalar r1, r2, z1, z2;
//
//    auto box = shape().boundingBox();
//    r1 = box.min_corner().x;
//    r2 = box.max_corner().x;
//    z1 = box.min_corner().y;
//    z2 = box.max_corner().y;
//
//    Scalar tstVol = M_PI * (r2*r2 - r1*r1) * (z2 - z1) / (2 * M_PI);
//
//    std::cout << volume / 3. << " = " << tstVol << "\n";

    return volume / 3.;
}

void Cell::addDiagonalLink(const Cell &cell)
{
    diagonalLinks_.push_back(CellLink(*this, cell));
    cellLinks_.push_back(CellLink(*this, cell));
}

void Cell::addBoundaryLink(const Face &face)
{
    if (!face.isBoundary())
        throw Exception("Cell", "addBoundaryLink", "cannot add a boundary link to a non-boundary face.");

    boundaryLinks_.push_back(BoundaryLink(*this, face));
}

void Cell::addInteriorLink(const Face &face, const Cell &cell)
{
    if (!face.isInterior())
        throw Exception("Cell", "addInteriorLink", "cannot add an interior link to a non-interior face.");

    interiorLinks_.push_back(InteriorLink(*this, face, cell));
    cellLinks_.push_back(CellLink(*this, cell));
}

const Cell &Cell::faceNeighbour(const Node &lNode, const Node &rNode) const
{
    for (const InteriorLink &nb: interiorLinks_)
    {
        const Face &face = nb.face();

        if ((lNode.id() == face.lNode().id() && rNode.id() == face.rNode().id()) ||
            (lNode.id() == face.rNode().id() && rNode.id() == face.lNode().id()))
        {
            return nb.cell();
        }
    }

    throw Exception("Cell", "faceNeighbour", "nodes do not belong to this cell.");
}

const std::vector<Ref<const Node> > Cell::nodes() const
{
    std:: vector<Ref<const Node>> nodes;
    nodes.reserve(nodeIds_.size());

    std::transform(nodeIds_.begin(), nodeIds_.end(), std::back_inserter(nodes), [this](Label id) {
        return std::cref(grid_.nodes()[id]);
    });

    return nodes;
}

bool Cell::isInCell(const Point2D &point) const
{
    return cellShape_.isInside(point);
}

//- Private methods

//- External functions
bool cellsShareFace(const Cell &cellA, const Cell &cellB)
{
    for (const InteriorLink &nb: cellA.neighbours())
    {
        if (nb.cell().id() == cellB.id())
            return true;
    }

    return false;
}
