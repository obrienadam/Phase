#ifndef FINITE_VOLUME_FIELD
#define FINITE_VOLUME_FIELD

#include <deque>

#include "Field.h"
#include "FiniteVolumeGrid2D.h"
#include "Input.h"
#include "SparseVector.h"

template<class T>
class FiniteVolumeField : public Field<T>
{
public:
    enum BoundaryType{FIXED, NORMAL_GRADIENT, SYMMETRY, OUTFLOW};

    FiniteVolumeField(const FiniteVolumeGrid2D& grid, const std::string& name);
    FiniteVolumeField(const Input& input, const FiniteVolumeGrid2D& grid, const std::string& name);
    FiniteVolumeField(const FiniteVolumeField& other);

    void fill(const T& val);
    void fillInterior(const T& val);

    void copyBoundaryTypes(const FiniteVolumeField& other);

    BoundaryType boundaryType(size_t faceId) const;
    T boundaryRefValue(size_t faceId) const;

    std::pair<BoundaryType, T> boundaryInfo(size_t faceId) const;

    const std::vector<T>& faces() const { return faces_; }
    std::vector<T>& faces() { return faces_; }

    FiniteVolumeField& savePreviousTimeStep(Scalar timeStep, int nPreviousFields);
    FiniteVolumeField& savePreviousIteration();

    FiniteVolumeField& prev(int i = 0) { return previousTimeSteps_[i].second; }
    const FiniteVolumeField& prev(int i = 0) const { return previousTimeSteps_[i].second; }
    Scalar prevTimeStep(int i = 0) { return previousTimeSteps_[i].first; }

    FiniteVolumeField& prevIter() { return previousIteration_.front(); }
    const FiniteVolumeField& prev() const { return previousIteration_.front(); }

    SparseVector sparseVector() const;

    FiniteVolumeField& operator=(const FiniteVolumeField& rhs);
    FiniteVolumeField& operator=(const SparseVector& rhs);
    FiniteVolumeField& operator+=(const FiniteVolumeField& rhs);
    FiniteVolumeField& operator-=(const FiniteVolumeField& rhs);
    FiniteVolumeField& operator*=(const FiniteVolumeField<Scalar>& rhs);
    FiniteVolumeField& operator*=(Scalar rhs);
    FiniteVolumeField& operator/=(Scalar lhs);
    FiniteVolumeField& operator/=(const FiniteVolumeField<Scalar>& rhs);

    const FiniteVolumeGrid2D& grid;

protected:

    void setBoundaryTypes(const Input& input);
    void setBoundaryRefValues(const Input& input);

    std::map<size_t, std::pair<BoundaryType, T> > patchBoundaries_;

    std::vector<T> faces_;

    std::deque< std::pair<Scalar, FiniteVolumeField<T> > > previousTimeSteps_;
    std::deque< FiniteVolumeField<T> > previousIteration_;
};

template<class T>
void interpolateFaces(FiniteVolumeField<T>& field);

template<class T>
void harmonicInterpolateFaces(FiniteVolumeField<T>& field);

template<class T>
FiniteVolumeField<T> smooth(const Field<T>& field, const std::vector< std::vector< Ref<const Cell> > >& rangeSearch, Scalar epsilon);

#include "FiniteVolumeField.tpp"

#endif
