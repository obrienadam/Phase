#ifndef PHASE_LINK_H
#define PHASE_LINK_H

#include "Geometry/Vector2D.h"

class Cell;

class Link
{
public:

    Link(const Cell &self);

    explicit Link(const Link &other);

    Link &operator=(const Link &rhs);

    const Cell &self() const
    { return self_; }

protected:

    const Cell &self_;
};

#endif
