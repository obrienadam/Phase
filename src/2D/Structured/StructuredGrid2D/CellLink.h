#ifndef PHASE_CELL_LINK_H
#define PHASE_CELL_LINK_H

#include "Geometry/Vector2D.h"

class Cell;

class CellLink {
public:
  CellLink() : _cellA(nullptr), _cellB(nullptr) {}

  CellLink(const Cell &cellA, const Cell &cellB);

  const Cell &cellA() const;

  const Cell &cellB() const;

  const Vector2D &r() const { return _rvec; }

protected:
  const Cell *_cellA, *_cellB;

  Vector2D _rvec;
};

#endif
