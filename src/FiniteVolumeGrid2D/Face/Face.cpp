#include <algorithm>

#include "Face.h"
#include "Exception.h"
#include "FiniteVolumeGrid2D.h"

Face::Face(Label lNodeId, Label rNodeId, const FiniteVolumeGrid2D &grid, Type type)
        :
        type_(type),
        nodes_(grid.nodes()),
        cells_(grid.cells()),
        nodeIds_(lNodeId, rNodeId)
{
    centroid_ = 0.5 * (lNode() + rNode());
    tangent_ = rNode() - lNode();
    normal_ = tangent_.normalVec();

    id_ = grid.faces().size();
}

Vector2D Face::polarOutwardNorm(const Point2D& point) const
{
    const Point2D &n1 = lNode();
    const Point2D &n2 = rNode();
    Vector2D sf((n2.y * n2.y - n1.y * n1.y) / 2., (n1.y + n2.y) / 2. * (n2.x - n1.x));
    return dot(sf, centroid_ - point) > 0. ? sf : -sf;
}

Vector2D Face::outwardNorm(const Point2D &point) const
{
    return dot(centroid_ - point, normal_) > 0. ? normal_ : -normal_;
}

Vector2D Face::outwardNorm() const
{
    return outwardNorm(lCell().centroid());
}

Scalar Face::volumeWeight() const
{
    Scalar v1 = rCell().volume();
    Scalar v2 = lCell().volume();
    return v1 / (v1 + v2);
}

Scalar Face::distanceWeight() const
{
    Scalar l1 = (centroid_ - rCell().centroid()).mag();
    Scalar l2 = (centroid_ - lCell().centroid()).mag();
    return l1 / (l1 + l2);
}

std::vector<Ref<const Cell> > Face::cells() const
{
    std::set<Label> ids;
    for(Label id: lNode().cellIds())
        ids.insert(id);

    for(Label id: rNode().cellIds())
        ids.insert(id);

    std::vector<Ref<const Cell>> cells;

    for(Label id: ids)
        cells.push_back(std::cref(cells_[id]));

    return cells;
}

void Face::addCell(const Cell &cell)
{
    if (type_ == INTERIOR)
    {
        if (cellIds_.size() == 2)
            throw Exception("Face", "addCell",
                            "an interior face cannot be shared between more than two cells. " + info());
    }
    else if (type_ == BOUNDARY)
    {
        if (cellIds_.size() == 1)
            throw Exception("Face", "addCell", "a boundary face cannot be shared by more than one cell. " + info());
    }

    cellIds_.push_back(cell.id());
}

std::string Face::info() const
{
    using namespace std;

    return "Face info: "
                   "Id: " + to_string(id_) + " "
                   "Node 1: " + to_string(lNode().id()) + ", "
                   "Node 2: " + to_string(rNode().id()) + ", "
                   "Centroid: " + centroid_.toString() + ", "
                   "Normal: " + normal_.toString();
}

//- Private methods
