#ifndef PATCH_H
#define PATCH_H

#include <string>
#include <vector>
#include <stdio.h>

#include "Types.h"

class Face;

class Patch
{
   public:

    Patch(const std::string& name, Label id) : name(name), id_(id) {}
    Patch(const Patch& other);

    size_t id() const { return id_; }

    void addFace(const Face& face);
    void addFaces(const std::vector< Ref<Face> >& faces);

    const std::vector< Ref<const Face> >& faces() const { return faces_; }

    std::string name;

private:

    size_t id_;

    std::vector< Ref<const Face> > faces_;

    friend class FiniteVolumeGrid2D;
};

#endif
