#include "Link.h"
#include "Cell.h"

Link::Link(const Cell &self)
    :
      self_(self)
{

}

Link::Link(const Link &other)
    :
      Link(other.self_)
{

}

const Cell& Link::self() const
{
    return self_;
}
