#ifndef VECTOR_FINITE_VOLUME_FIELD_H
#define VECTOR_FINITE_VOLUME_FIELD_H

#include "Field.h"
#include "Vector2D.h"
#include "SparseVector.h"
#include "Input.h"

class FiniteVolumeGrid2D;

class VectorFiniteVolumeField : public Field<Vector2D>
{
public:

    enum BoundaryType{FIXED, NORMAL_GRADIENT};

    VectorFiniteVolumeField(const FiniteVolumeGrid2D& grid, const std::string &name);
    VectorFiniteVolumeField(const Input& input, const FiniteVolumeGrid2D& grid, const std::string& name);

    VectorFiniteVolumeField& operator =(const SparseVector& rhs);

    std::vector<Vector2D>& faces() { return faces_; }
    const std::vector<Vector2D>& faces() const { return faces_; }

    BoundaryType boundaryType(size_t faceId) const;

    const FiniteVolumeGrid2D &grid;

protected:

    std::vector<BoundaryType> boundaryTypes_;
    std::vector<Vector2D> boundaryRefValues_;

    std::vector<Vector2D> faces_;

};

#endif
