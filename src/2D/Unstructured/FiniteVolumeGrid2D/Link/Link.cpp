#include "FiniteVolumeGrid2D/Cell/Cell.h"

#include "Link.h"

Link::Link(const Cell &self) : self_(self) {}

Link::Link(const Link &other) : Link(other.self_) {}
