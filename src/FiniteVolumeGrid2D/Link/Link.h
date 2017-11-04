#ifndef LINK_H
#define LINK_H

#include "Vector2D.h"
#include "Types.h"

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
