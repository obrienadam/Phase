#include <algorithm>

#include "Face.h"
#include "Cell.h"
#include "Exception.h"
#include "FiniteVolumeGrid2D.h"

Face::Face(Label lNodeId, Label rNodeId, const FiniteVolumeGrid2D &grid, Type type)
    :
      type_(type),
      patchPtr_(nullptr),
      nodes_(grid.nodes()),
      cells_(grid.cells()),
      nodeIds_(lNodeId, rNodeId)
{
    centroid_ = 0.5*(lNode() + rNode());
    tangent_ = rNode() - lNode();
    normal_ = tangent_.normalVec();

    id_ = grid.faces().size();
}

Vector2D Face::outwardNorm(const Point2D& point) const
{
    return dot(centroid_ - point, normal_) > 0. ? normal_ : -normal_;
}

void Face::addCell(const Cell &cell)
{
    if(type_ == INTERIOR)
    {
        if(cellIds_.size() == 2)
            throw Exception("Face", "addCell", "an interior face cannot be shared between more than two cells. " + info());
    }
    else if(type_ == BOUNDARY)
    {
            if(cellIds_.size() == 1)
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

void Face::addToPatch(const Patch& patch) const
{
    if(patchPtr_ != nullptr)
        throw Exception("Face", "addToPatch", "face already belongs to a patch.");

    patchPtr_ = &patch;
}

void Face::changePatch(const Patch &patch) const
{
    patchPtr_ = &patch;
}
