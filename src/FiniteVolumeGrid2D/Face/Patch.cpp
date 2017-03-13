#include "Patch.h"
#include "FiniteVolumeGrid2D.h"

Patch::Patch(const Patch &other)
        :
        Patch(other.name, other.id_)
{
    faces_ = other.faces_;

    for (const Face &face: faces_)
        face.changePatch(*this);
}

void Patch::addFace(const Face &face)
{
    if (face.isBoundary())
    {
        faces_.push_back(Ref<const Face>(face));
        face.addToPatch(*this);
    }
}

void Patch::addFaces(const std::vector<Ref<Face> > &faces)
{
    for (Face &face: faces)
    {
        if (face.isBoundary())
        {
            faces_.push_back(Ref<const Face>(face));
            face.addToPatch(*this);
        }
    }
}
