#ifndef PHASE_SCALAR_GRADIENT_H
#define PHASE_SCALAR_GRADIENT_H

#include "VectorFiniteVolumeField.h"

class ScalarGradient : public VectorFiniteVolumeField {
public:
  enum Method { FACE_TO_CELL, GREEN_GAUSS_CELL, GREEN_GAUSS_NODE };

  static Vector2D computeGradient(const ScalarFiniteVolumeField &phi,
                                  const Cell &c);

  explicit ScalarGradient(const ScalarFiniteVolumeField &phi,
                          const std::shared_ptr<const CellGroup> &cells);

  void computeFaces();

  void compute(const CellGroup &cells, Method method = FACE_TO_CELL);

  void compute(Method method = FACE_TO_CELL);

  void computeAxisymmetric(const CellGroup &cells);

  void computeAxisymmetric(const ScalarFiniteVolumeField &cw,
                           const ScalarFiniteVolumeField &fw,
                           const CellGroup &cells);

private:
  const ScalarFiniteVolumeField &phi_;
};

#endif
