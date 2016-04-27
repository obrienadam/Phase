#ifndef LINK_H
#define LINK_H

#include "Vector2D.h"
#include "Types.h"

class Cell;
class Face;

class Link
{
public:

    Link(const Cell& self) : self_(self) {}

    const Cell& self() const { return self_; }

protected:

    Ref<const Cell> self_;
};

#endif
