#ifndef PHASE_FINITE_VOLUME_FIELD
#define PHASE_FINITE_VOLUME_FIELD

#include "System/Input.h"

#include "Field.h"
#include "FiniteVolumeGrid2D/FiniteVolumeGrid2D.h"
#include "FiniteVolume/Equation/IndexMap.h"

template<class T>
class FiniteVolumeField : public Field<T>
{
public:
    enum BoundaryType
    {
        FIXED, NORMAL_GRADIENT, SYMMETRY, OUTFLOW, PARTIAL_SLIP
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
                               const std::shared_ptr<const CellGroup> &cellGroup = nullptr,
                               const std::shared_ptr<IndexMap> &indexMap = nullptr);

    explicit FiniteVolumeField(const Input &input,
                               const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                               const std::string &name,
                               const T &val = T(),
                               bool faces = true,
                               bool nodes = false,
                               const std::shared_ptr<const CellGroup> &cellGroup = nullptr,
                               const std::shared_ptr<IndexMap> &indexMap = nullptr);

    //- Initialization
    void fill(const T &val);

    void fill(const T &val, const CellGroup &group);

    void fillInterior(const T &val);

    void assign(const FiniteVolumeField<T> &field);

    std::shared_ptr<IndexMap> &indexMap()
    { return indexMap_; }

    const std::shared_ptr<IndexMap> &indexMap() const
    { return indexMap_; }

    template<class TFunc>
    void computeCells(const TFunc &fcn)
    {
        for (const Cell &cell: grid_->cells())
            (*this)(cell) = fcn(cell);
    }

    template<class TFunc>
    void computeFaces(const TFunc &fcn)
    {
        for (const Face &face: grid_->faces())
            (*this)(face) = fcn(face);
    }

    template<class TFunc>
    void computeInteriorFaces(const TFunc &fcn)
    {
        for (const Face &face: grid_->interiorFaces())
            (*this)(face) = fcn(face);
    }

    template<class TFunc>
    void computeBoundaryFaces(const TFunc &fcn)
    {
        for (const Face &face: grid_->boundaryFaces())
            (*this)(face) = fcn(face);
    }

    void faceToCell(const CellGroup &cells);

    template<class UnaryPredicate>
    void faceToCell(const FiniteVolumeField<Scalar> &cellWeight,
                    const FiniteVolumeField<Scalar> &faceWeight,
                    const CellGroup &cells,
                    const UnaryPredicate &p);

    void faceToCell(const FiniteVolumeField<Scalar> &cellWeight,
                    const FiniteVolumeField<Scalar> &faceWeight,
                    const CellGroup &cells);

    void faceToCellAxisymmetric(const CellGroup &cells);

    void faceToCellAxisymmetric(const FiniteVolumeField<Scalar> &cw, const FiniteVolumeField<Scalar> &fw, const CellGroup &cells);

    //- Boundaries
    void copyBoundaryTypes(const FiniteVolumeField &other);

    BoundaryType boundaryType(const FaceGroup &patch) const;

    BoundaryType boundaryType(const Face &face) const;

    T boundaryRefValue(const FaceGroup &patch) const;

    T boundaryRefValue(const Face &face) const;

    template<class TFunc>
    void interpolateFaces(const TFunc &alpha);

    void interpolateFaces(InterpolationType type = DISTANCE);

    void setBoundaryFaces();

    void setBoundaryFaces(BoundaryType bType, const std::function<T(const Face &face)> &fcn);

    void interpolateNodes();

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
    { return cellGroup_ ? *cellGroup_ : grid_->localCells(); }

    virtual void setCellGroup(const std::shared_ptr<const CellGroup> &cellGroup)
    { cellGroup_ = cellGroup; }

    virtual void setIndexMap(const std::shared_ptr<IndexMap> &indexMap)
    { indexMap_ = indexMap; }

    //- Field history
    FiniteVolumeField &savePreviousTimeStep(Scalar timeStep, int nPreviousFields);

    FiniteVolumeField &savePreviousIteration();

    void clearHistory();

    FiniteVolumeField &oldField(int i)
    { return *previousTimeSteps_[i].second; }

    const FiniteVolumeField &oldField(int i) const
    { return *previousTimeSteps_[i].second; }

    Scalar oldTimeStep(int i) const
    { return previousTimeSteps_[i].first; }

    const FiniteVolumeField &prevIteration() const
    { return *previousIteration_; }

    //- Parallel

    void sendMessages();

    //- Operators

    FiniteVolumeField &operator+=(const FiniteVolumeField &rhs);

    FiniteVolumeField &operator-=(const FiniteVolumeField &rhs);

    FiniteVolumeField &operator*=(const FiniteVolumeField<Scalar> &rhs);

    FiniteVolumeField &operator/=(const FiniteVolumeField<Scalar> &rhs);

    FiniteVolumeField &operator*=(Scalar rhs);

    FiniteVolumeField &operator/=(Scalar lhs);

    void setGrid(const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    const std::shared_ptr<const FiniteVolumeGrid2D> &grid() const
    { return grid_; }

    //- Debug
    void writeToFile(const std::string &filename) const;

    bool isfinite() const;

protected:

    typedef std::pair<Scalar, FiniteVolumeField<T>> PreviousField;

    //- Parallel
    static std::vector<std::vector<T>> sendBuffers_, recvBuffers_;

    void setBoundaryTypes(const Input &input);

    void setBoundaryRefValues(const Input &input);

    //- Data members
    std::unordered_map<std::string, std::pair<BoundaryType, T> > patchBoundaries_;

    //- Grid
    std::shared_ptr<const FiniteVolumeGrid2D> grid_;

    //- Main cell group
    std::shared_ptr<const CellGroup> cellGroup_;

    //- Misc data
    std::vector<T> faces_, nodes_;

    //- Field history
    std::deque<std::pair<Scalar, std::shared_ptr<FiniteVolumeField<T>>>> previousTimeSteps_;

    std::shared_ptr<FiniteVolumeField<T>> previousIteration_;

    //- Index map
    std::shared_ptr<IndexMap> indexMap_;
};

#include "FiniteVolumeField.tpp"

#endif
