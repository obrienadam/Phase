#include "Cell.h"
#include "Polygon.h"
#include "Exception.h"

Cell::Cell(const std::vector<size_t> &faceIds, std::vector<Face> &faces, bool isActive)
{   
    isActive_ = isActive;

    for (size_t faceId: faceIds)
    {
        Ref<const Face> face = faces[faceId];
        faces_.push_back(face);
    }

    std::vector<Point2D> vertices;

    for(int i = 0, end = faces_.size(); i < end; ++i)
    {
        const Face &face1 = faces_[i], &face2 = faces_[(i + 1)%faces_.size()];

        if (face1.rNode() == face2.lNode() || face1.rNode() == face2.rNode())
        {
            vertices.push_back(face1.rNode());
            nodeIds_.push_back(face1.rNode().id());
        }
        else if(face1.lNode() == face2.lNode() || face1.lNode() == face2.rNode())
        {
            vertices.push_back(face1.lNode());
            nodeIds_.push_back(face1.lNode().id());
        }
        else
        {
            throw Exception("Cell", "Cell", "adjacent faces do not share a node.");
        }
    }

    cellShape_ = Polygon(vertices);

    if(!cellShape_.isConvex())
        throw Exception("Cell", "Cell", "non-convex cells are not allowed.");

    volume_ = cellShape_.area();
    centroid_ = cellShape_.centroid();

    if(volume_ < 0.)
        throw Exception("Cell", "Cell", "faces are not oriented in a counter-clockwise manner.");
}

bool Cell::isInCell(const Point2D &point) const
{
    return cellShape_.isInside(point);
}

//- Private methods

void Cell::computeCellAdjacency()
{
    interiorLinks_.clear();
    boundaryLinks_.clear();

    for(const Face& face: faces_)
    {
        if(face.isInterior())
        {
            const Cell& cell = this == &face.lCell() ? face.rCell() : face.lCell();
            InteriorLink link(*this, face, cell);
            interiorLinks_.push_back(link);
        }
        else if(face.isBoundary())
        {
            BoundaryLink link(*this, face);
            boundaryLinks_.push_back(link);
        }
        else
        {
            throw Exception("Cell", "computeCellAdjacency", "unrecognized face type.");
        }
    }
}
