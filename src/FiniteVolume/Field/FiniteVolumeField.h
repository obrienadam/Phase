#ifndef FINITE_VOLUME_FIELD
#define FINITE_VOLUME_FIELD

#include "Field.h"
#include "FiniteVolumeGrid2D.h"
#include "Input.h"
#include "Vector.h"

template<class T>
class FiniteVolumeField : public Field<T>
{
public:
    enum BoundaryType
    {
        FIXED, NORMAL_GRADIENT, SYMMETRY, OUTFLOW
    };

    enum InterpolationType
    {
        VOLUME, DISTANCE
    };

    //- Constructors
    explicit FiniteVolumeField(const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                               const std::string &name,
                               const T &val = T(),
                               bool faces = true,
                               bool nodes = false,
                               const std::shared_ptr<const CellGroup> &cellGroup = nullptr);

    explicit FiniteVolumeField(const Input &input,
                               const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                               const std::string &name,
                               const T &val = T(),
                               bool faces = true,
                               bool nodes = false,
                               const std::shared_ptr<const CellGroup> &cellGroup = nullptr);

    //- Initialization
    void fill(const T &val);

    void fillInterior(const T &val);

    template<class TFunc>
    void computeCells(const TFunc &fcn)
    {
        for (const Cell &cell: grid().cells())
            (*this)(cell) = fcn(cell);
    }

    template<class TFunc>
    void computeFaces(const TFunc &fcn)
    {
        for (const Face &face: grid().faces())
            (*this)(face) = fcn(face);
    }

    template<class TFunc>
    void computeInteriorFaces(const TFunc &fcn)
    {
        for (const Face &face: grid().interiorFaces())
            (*this)(face) = fcn(face);
    }

    template<class TFunc>
    void computeBoundaryFaces(const TFunc &fcn)
    {
        for (const Face &face: grid().boundaryFaces())
            (*this)(face) = fcn(face);
    }

    void faceToCell(const FiniteVolumeField<Scalar> &cellWeight,
                    const FiniteVolumeField<Scalar> &faceWeight,
                    const CellGroup &cells);

    void setPatch(const Patch &patch, const std::function<T(const Face &face)> &fcn)
    {
        for (const Face &face: patch)
            (*this)(face) = fcn(face);
    }

    //- Boundaries
    void copyBoundaryTypes(const FiniteVolumeField &other);

    BoundaryType boundaryType(const Patch &patch) const;

    BoundaryType boundaryType(const Face &face) const;

    T boundaryRefValue(const Patch &patch) const;

    std::pair<BoundaryType, T> boundaryInfo(const Face &face) const;

    template<class TFunc>
    void interpolateFaces(const TFunc &alpha)
    {
        auto &self = *this;

        for (const Face &face: grid_->interiorFaces())
        {
            Scalar g = alpha(face);
            self(face) = g * self(face.lCell()) + (1. - g) * self(face.rCell());
        }

        setBoundaryFaces();
    }

    void interpolateFaces(InterpolationType type = VOLUME)
    {
        switch (type)
        {
            case VOLUME:
                interpolateFaces([](const Face &face) {
                    return face.volumeWeight();
                });
                break;
            case DISTANCE:
                interpolateFaces([](const Face &face) {
                    return face.distanceWeight();
                });
                break;
        }
    }

    void setBoundaryFaces();

    void setBoundaryFaces(BoundaryType bType, const std::function<T(const Face &face)> &fcn);

    //- Field info
    bool hasFaces() const
    { return !faces_.empty(); }

    bool hasNodes() const
    { return !nodes_.empty(); }

    //- Face-centered values
    const std::vector<T> &faces() const
    { return faces_; }

    std::vector<T> &faces()
    { return faces_; }

    //- Node-centered values
    std::vector<T> &nodes()
    { return nodes_; }

    const std::vector<T> &nodes() const
    { return nodes_; }

    //- Access operators
    T &operator()(const Cell &cell)
    { return std::vector<T>::operator[](cell.id()); }

    const T &operator()(const Cell &cell) const
    { return std::vector<T>::operator[](cell.id()); }

    T &operator()(Label id)
    { return std::vector<T>::operator[](id); }

    const T &operator()(Label id) const
    { return std::vector<T>::operator[](id); }

    T &operator()(const Face &face)
    { return faces_[face.id()]; }

    const T &operator()(const Face &face) const
    { return faces_[face.id()]; }

    T &operator()(const Node &node)
    { return nodes_[node.id()]; }

    const T &operator()(const Node &node) const
    { return nodes_[node.id()]; }

    //- Cell group access (by default returns local active cells)
    const CellGroup &cells() const
    { return cellGroup_ ? *cellGroup_ : grid_->localActiveCells(); }

    void setCellGroup(const CellGroup& cellGroup)
    { cellGroup_ = std::make_shared<CellGroup>(cellGroup); }

    void setCellGroup(CellGroup&& cellGroup)
    { cellGroup_ = std::make_shared<CellGroup>(cellGroup); }

    void setCellGroup(const std::shared_ptr<const CellGroup>& cellGroup)
    { cellGroup_ = cellGroup; }

    //- Field history
    FiniteVolumeField &savePreviousTimeStep(Scalar timeStep, int nPreviousFields);

    FiniteVolumeField &savePreviousIteration();

    void clearHistory();

    FiniteVolumeField &oldField(int i)
    { return previousTimeSteps_[i]->second; }

    const FiniteVolumeField &oldField(int i) const
    { return previousTimeSteps_[i]->second; }

    Scalar oldTimeStep(int i) const
    { return previousTimeSteps_[i]->first; }

    const FiniteVolumeField &prevIteration() const
    { return *previousIteration_; }

    //- Vectorization
    Vector vectorize() const;

    //- Operators
    FiniteVolumeField &operator=(const Vector &rhs);

    FiniteVolumeField &operator+=(const FiniteVolumeField &rhs);

    FiniteVolumeField &operator-=(const FiniteVolumeField &rhs);

    FiniteVolumeField &operator*=(const FiniteVolumeField<Scalar> &rhs);

    FiniteVolumeField &operator/=(const FiniteVolumeField<Scalar> &rhs);

    FiniteVolumeField &operator*=(Scalar rhs);

    FiniteVolumeField &operator/=(Scalar lhs);

    std::shared_ptr<const FiniteVolumeGrid2D> gridPtr() const
    { return grid_; }

    const FiniteVolumeGrid2D &grid() const
    { return *grid_; }

    //- Debug
    void writeToFile(const std::string &filename) const;

protected:

    typedef std::pair<Scalar, FiniteVolumeField<T>> PreviousField;

    void setBoundaryTypes(const Input &input);

    void setBoundaryRefValues(const Input &input);

    std::map<Label, std::pair<BoundaryType, T> > patchBoundaries_;

    //- Grid
    std::shared_ptr<const FiniteVolumeGrid2D> grid_;

    //- Main cell group
    std::shared_ptr<const CellGroup> cellGroup_;

    //- Misc data
    std::vector<T> faces_, nodes_;

    //- Field history
    std::vector<std::shared_ptr<PreviousField>> previousTimeSteps_;

    std::shared_ptr<FiniteVolumeField<T>> previousIteration_;
};

#include "FiniteVolumeField.tpp"

#endif
