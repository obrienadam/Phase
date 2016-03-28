#ifndef PATCH_H
#define PATCH_H

#include <string>
#include <vector>
#include <stdio.h>

#include "Types.h"
#include "BoundingBox.h"

class FiniteVolumeGrid2D;
class Face;

class Patch
{
   public:

    Patch(size_t id, const std::string& name = "Unnamed_patch") : id_(id), name(name) {}
    Patch(const Patch& other);

    size_t id() const { return id_; }

    void addFace(Face& face);
    void addFaces(FiniteVolumeGrid2D& grid, const BoundingBox& bBox);
    void addFaces(const std::vector< Ref<Face> >& faces);

    const std::vector< Ref<const Face> >& faces() const { return faces_; }

    std::string name;

private:

    size_t id_;

    std::vector< Ref<const Face> > faces_;

    friend class FiniteVolumeGrid2D;
};

#endif
