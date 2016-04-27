#include <algorithm>

#include "Face.h"
#include "Cell.h"
#include "Exception.h"

Face::Face(size_t lNodeId, size_t rNodeId, const std::vector<Node> &nodes, Type type)
    :
      type_(type),
      patchPtr_(nullptr)
{
    nodes_.push_back(nodes[lNodeId]);
    nodes_.push_back(nodes[rNodeId]);

    centroid_ = 0.5*(nodes_[0] + nodes_[1]);
    tangent_ = nodes_[1] - nodes_[0];
    normal_ = tangent_.normalVec();
}

Vector2D Face::outwardNorm(const Point2D& point) const
{
    return dot(centroid_ - point, normal_) > 0. ? normal_ : -normal_;
}

void Face::addCell(const Cell &cell)
{
    if(type_ == INTERIOR)
    {
        if(cells_.size() == 2)
            throw Exception("Face", "addCell", "an interior face cannot be shared between more than two cells. " + info());
    }
    else if(type_ == BOUNDARY)
    {
            if(cells_.size() == 1)
                throw Exception("Face", "addCell", "a boundary face cannot be shared by more than one cell. " + info());
    }

    cells_.push_back(cell);
}

std::string Face::info() const
{
    using namespace std;

    return "Face info: "
            "Id: " + to_string(id_) + " "
            "Centroid: " + centroid_.toString() + " "
            "Normal: " + normal_.toString() + " "
            "Tangent: " + normal_.toString() + " ";
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
