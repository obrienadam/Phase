#ifndef PHASE_CELL_H
#define PHASE_CELL_H

#include "BoundaryFaceStencil.h"
#include "Geometry/RectangularPrism.h"
#include "InteriorFaceStencil.h"

class StructuredGrid3D;

class Cell {
public:
  enum Index { I, J, K };

  Cell(const StructuredGrid3D &grid, Label i, Label j, Label k);

  const StructuredGrid3D &grid() const { return _grid; }

  Label i() const { return _i; }

  Label j() const { return _j; }

  Label k() const { return _k; }

  Label id() const { return _id; }

  const Cell &nb(Index idx, int offset) const;

  Scalar volume() const { return _shape.volume(); }

  const Point3D &centroid() const { return _shape.centroid(); }

  //- Get stencils
  const std::vector<InteriorFaceStencil> &interiorStencils() const {
    return _interiorStencils;
  }

  const std::vector<BoundaryFaceStencil> &boundaryStencils() const {
    return _boundaryStencils;
  }

  //- Init stencils
  void initStencils(int order);

protected:
  std::vector<InteriorFaceStencil> _interiorStencils;

  std::vector<BoundaryFaceStencil> _boundaryStencils;

  const StructuredGrid3D &_grid;

  Label _i, _j, _k, _id;

  RectangularPrism _shape;
};

#endif
