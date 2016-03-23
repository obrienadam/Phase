#include "Patch.h"
#include "FiniteVolumeGrid2D.h"

Patch::Patch(const Patch &other)
    :
      Patch(other.id_, other.name)
{
    faces_ = other.faces_;

    for(const Face &face: faces_)
        face.changePatch(*this);
}

void Patch::addFace(Face &face)
{
    if(face.isBoundary())
    {
        faces_.push_back(Ref<const Face>(face));
        face.addToPatch(*this);
    }
}

void Patch::addFaces(FiniteVolumeGrid2D &grid, const BoundingBox &bBox)
{
    for(Face& face: grid.faces)
    {
        if(face.isBoundary() && bBox.isInBox(face.centroid()))
        {
            faces_.push_back(Ref<const Face>(face));
            face.addToPatch(*this);
        }
    }
}

void Patch::addFaces(const std::vector<Ref<Face> > &faces)
{
    for(Face& face: faces)
    {
        if(face.isBoundary())
        {
            faces_.push_back(Ref<const Face>(face));
            face.addToPatch(*this);
        }
    }
}
