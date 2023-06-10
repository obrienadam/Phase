#include "BoundaryPatch.h"

BoundaryPatch::BoundaryPatch(const std::string &name) : Set<Face>(name) {}

bool BoundaryPatch::add(const Face &face) {
  if (face.isBoundaryFace())
    Set<Face>::add(face);
}
