#ifndef SCALAR_FINITE_VOLUME_FIELD_H
#define SCALAR_FINITE_VOLUME_FIELD_H

#include <string>

#include "Field.h"
#include "SparseVector.h"
#include "Input.h"

class FiniteVolumeGrid2D;

class ScalarFiniteVolumeField : public Field<Scalar>
{
public:

    enum BoundaryType{FIXED, NORMAL_GRADIENT};

    ScalarFiniteVolumeField(const FiniteVolumeGrid2D& grid, const std::string& name);
    ScalarFiniteVolumeField(const Input& input, const FiniteVolumeGrid2D& grid, const std::string& name);

    void fill(Scalar val);

    ScalarFiniteVolumeField& operator =(const SparseVector& rhs);

    const std::vector<Scalar>& faces() const { return faces_; }
    std::vector<Scalar>& faces() { return faces_; }

    BoundaryType boundaryType(size_t faceId) const;

    ScalarFiniteVolumeField& operator*=(const ScalarFiniteVolumeField& rhs);

    const FiniteVolumeGrid2D& grid;

protected:

    std::vector<BoundaryType> boundaryTypes_;
    std::vector<Scalar> boundaryRefValues_;

    std::vector<Scalar> faces_;
};

ScalarFiniteVolumeField operator*(const ScalarFiniteVolumeField& lhs, ScalarFiniteVolumeField rhs);

#endif
