#ifndef FINITE_VOLUME_FIELD
#define FINITE_VOLUME_FIELD

#include <memory>

#include "Field.h"
#include "FiniteVolumeGrid2D.h"
#include "Input.h"
#include "SparseVector.h"

template<class T>
class FiniteVolumeField : public Field<T>
{
public:
    enum BoundaryType{FIXED, NORMAL_GRADIENT};

    FiniteVolumeField(const FiniteVolumeGrid2D& grid, const std::string& name);
    FiniteVolumeField(const Input& input, const FiniteVolumeGrid2D& grid, const std::string& name);
    FiniteVolumeField(const FiniteVolumeField& other);

    void fill(const T& val);
    void fillInterior(const T& val);

    void copyBoundaryTypes(const FiniteVolumeField& other);

    BoundaryType boundaryType(size_t faceId) const;
    T boundaryRefValue(size_t faceId) const;

    const std::vector<T>& faces() const { return faces_; }
    std::vector<T>& faces() { return faces_; }

    FiniteVolumeField& save();
    FiniteVolumeField& prev() { return *prevFieldPtr_; }
    const FiniteVolumeField& prev() const { return *prevFieldPtr_; }

    FiniteVolumeField& operator=(const FiniteVolumeField& rhs);
    FiniteVolumeField& operator=(const SparseVector& rhs);
    FiniteVolumeField& operator+=(const FiniteVolumeField& rhs);
    FiniteVolumeField& operator-=(const FiniteVolumeField& rhs);
    FiniteVolumeField& operator*=(const FiniteVolumeField<Scalar>& rhs);
    FiniteVolumeField& operator*=(Scalar rhs);
    FiniteVolumeField& operator/=(Scalar lhs);

    const FiniteVolumeGrid2D& grid;

protected:

    void setBoundaryTypes(const Input& input);
    void setBoundaryRefValues(const Input& input);

    std::vector<BoundaryType> boundaryTypes_;
    std::vector<T> boundaryRefValues_;

    std::vector<T> faces_;

    std::shared_ptr<FiniteVolumeField> prevFieldPtr_;
};

template<class T>
void interpolateFaces(FiniteVolumeField<T>& field);

#include "FiniteVolumeField.tpp"

#endif
