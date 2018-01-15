#include "FiniteVolumeField.h"
#include "Exception.h"

#include <boost/algorithm/string.hpp>
#include <fstream>

//- Constructors

template<class T>
FiniteVolumeField<T>::FiniteVolumeField(const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                        const std::string &name,
                                        const T &val,
                                        bool faces,
                                        bool nodes,
                                        const std::shared_ptr<const CellGroup>& cellGroup)
        :
        Field<T>::Field(grid->cells().size(), val, name),
        grid_(grid)
{
    cellGroup_ = cellGroup;

    if (faces)
        faces_.resize(grid_->faces().size(), val);

    if (nodes)
        nodes_.resize(grid_->nodes().size(), val);
}

template<class T>
FiniteVolumeField<T>::FiniteVolumeField(const Input &input,
                                        const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                        const std::string &name,
                                        const T &val,
                                        bool faces,
                                        bool nodes,
                                        const std::shared_ptr<const CellGroup>& cellGroup)
        :
        FiniteVolumeField(grid, name, val, faces, nodes, cellGroup)
{
    setBoundaryTypes(input);
    setBoundaryRefValues(input);
}

//- Public methods

template<class T>
void FiniteVolumeField<T>::fill(const T &val)
{
    std::fill(this->begin(), this->end(), val);
    std::fill(faces_.begin(), faces_.end(), val);
    std::fill(nodes_.begin(), nodes_.end(), val);
}

template<class T>
void FiniteVolumeField<T>::fillInterior(const T &val)
{
    std::fill(this->begin(), this->end(), val);

    if (!faces_.empty())
    {
        for (const Face &face: grid().interiorFaces())
            faces_[face.id()] = val;
    }

    if (!nodes_.empty())
    {
        for (const Node &node: grid().interiorNodes())
            nodes_[node.id()] = val;
    }
}

template<class T>
void FiniteVolumeField<T>::assign(const FiniteVolumeField<T>& field)
{
    this->assign(field.begin(), field.end());
    faces_.assign(field.faces_.begin(), field.faces_.end());
    nodes_.assign(field.nodes_.begin(), field.nodes_.end());
    patchBoundaries_ = field.patchBoundaries_;
    grid_ = field.grid_;
    cellGroup_ = field.cellGroup_;
}

template<class T>
void FiniteVolumeField<T>::copyBoundaryTypes(const FiniteVolumeField &other)
{
    patchBoundaries_ = other.patchBoundaries_;
}

template<class T>
typename FiniteVolumeField<T>::BoundaryType FiniteVolumeField<T>::boundaryType(const Patch &patch) const
{
    auto it = patchBoundaries_.find(patch.id());
    return it == patchBoundaries_.end() ? NORMAL_GRADIENT : it->second.first;
}

template<class T>
typename FiniteVolumeField<T>::BoundaryType FiniteVolumeField<T>::boundaryType(const Face &face) const
{
    return boundaryType(grid().patch(face));
}

template<class T>
T FiniteVolumeField<T>::boundaryRefValue(const Patch &patch) const
{
    return patchBoundaries_.find(patch.id())->second.second;
}

template<class T>
void FiniteVolumeField<T>::setBoundaryFaces()
{
    auto &self = *this;

    for (const Patch &patch: grid().patches())
    {
        switch (boundaryType(patch))
        {
            case FIXED:
                break;
            case NORMAL_GRADIENT:
            case SYMMETRY:
                for(const Face& face: patch)
                    faces_[face.id()] = self[face.lCell().id()];
                break;
        }
    }
}

template<class T>
void FiniteVolumeField<T>::setBoundaryFaces(BoundaryType bType, const std::function<T(const Face &face)> &fcn)
{
    auto &self = *this;
    for (const Patch &patch: grid_->patches())
    {
        if (boundaryType(patch) == bType)
            for (const Face &face: patch)
                self(face) = fcn(face);
    }
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::savePreviousTimeStep(Scalar timeStep, int nPreviousFields)
{
    if(previousTimeSteps_.size() == nPreviousFields)
    {
        auto prevTimeStep = previousTimeSteps_.back();
        prevTimeStep->second = *this;
        prevTimeStep->second.clearHistory();
        previousTimeSteps_.insert(previousTimeSteps_.begin(), prevTimeStep);
        previousTimeSteps_.pop_back();
    }
    else
    {
        auto prevTimeStep = std::make_shared<PreviousField>(timeStep, *this);
        prevTimeStep->second.clearHistory();
        previousTimeSteps_.insert(previousTimeSteps_.begin(), prevTimeStep);
    }

    return previousTimeSteps_.front()->second;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::savePreviousIteration()
{
    previousIteration_ = std::make_shared<FiniteVolumeField<T>>(*this);
    previousIteration_->clearHistory();

    return *previousIteration_;
}

template<class T>
void FiniteVolumeField<T>::clearHistory()
{
    previousIteration_ = nullptr;
    previousTimeSteps_.clear();
}

//- Operators

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator+=(const FiniteVolumeField &rhs)
{
    auto &self = *this;

    for (const Cell &cell: grid().cells())
        self(cell) += rhs(cell);

    if (!faces_.empty() && rhs.hasFaces())
        for (const Face &face: grid().faces())
            self(face) += rhs(face);

    if (!nodes_.empty() && rhs.hasNodes())
        for (const Node &node: grid().nodes())
            self(node) += rhs(node);

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator-=(const FiniteVolumeField &rhs)
{
    auto &self = *this;

    for (const Cell &cell: grid().cells())
        self(cell) -= rhs(cell);

    if (!faces_.empty() && rhs.hasFaces())
        for (const Face &face: grid().faces())
            self(face) -= rhs(face);

    if (!nodes_.empty() && rhs.hasNodes())
        for (const Node &node: grid().nodes())
            self(node) -= rhs(node);

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator*=(const FiniteVolumeField<Scalar> &rhs)
{
    auto &self = *this;

    for (const Cell &cell: grid().cells())
        self(cell) *= rhs(cell);

    if (!faces_.empty() && rhs.hasFaces())
        for (const Face &face: grid().faces())
            self(face) *= rhs(face);

    if (!nodes_.empty() && rhs.hasNodes())
        for (const Node &node: grid().nodes())
            self(node) *= rhs(node);

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator/=(const FiniteVolumeField<Scalar> &rhs)
{
    auto &self = *this;

    for (const Cell &cell: grid().cells())
        self(cell) /= rhs(cell);

    if (hasFaces() && rhs.hasFaces())
        for (const Face &face: grid().faces())
            self(face) /= rhs(face);

    if (hasNodes() && rhs.hasNodes())
        for (const Node &node: grid().nodes())
            self(node) /= rhs(node);

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator*=(Scalar rhs)
{
    auto &self = *this;

    for (const Cell &cell: grid().cells())
        self(cell) *= rhs;

    if (!faces_.empty())
        for (const Face &face: grid().faces())
            self(face) *= rhs;

    if (!nodes_.empty())
        for (const Node &node: grid().nodes())
            self(node) *= rhs;

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator/=(Scalar rhs)
{
    auto &self = *this;

    for (const Cell &cell: grid().cells())
        self(cell) /= rhs;

    if (!faces_.empty())
        for (const Face &face: grid().faces())
            self(face) /= rhs;

    if (!nodes_.empty())
        for (const Node &node: grid().nodes())
            self(node) /= rhs;

    return self;
}

//- Debug

template<class T>
void FiniteVolumeField<T>::writeToFile(const std::string &filename) const
{
    std::ofstream fout(filename);

    for (const T &val: *this)
        fout << val << "\n";

    fout.close();
}

//- Protected methods

template<class T>
void FiniteVolumeField<T>::setBoundaryTypes(const Input &input)
{
    using namespace std;


    std::string typeStr = input.boundaryInput().get<string>("Boundaries." + Field<T>::name() + ".*.type", "");
    BoundaryType boundaryType;
    std::vector<Scalar> coeffs;

    if (!typeStr.empty())
    {
        if (typeStr == "fixed")
            boundaryType = FIXED;
        else if (typeStr == "normal_gradient")
            boundaryType = NORMAL_GRADIENT;
        else if (typeStr == "symmetry")
            boundaryType = SYMMETRY;
        else if (typeStr == "outflow")
            boundaryType = OUTFLOW;
        else
            throw Exception("FiniteVolumeField<T>", "setBoundaryTypes", "invalid boundary type \"" + typeStr + "\".");

        for (const Patch &patch: grid().patches())
        {
            patchBoundaries_[patch.id()] = std::make_pair(boundaryType, T());
        }
    }

    for (const Patch &patch: grid().patches())
    {
        typeStr = input.boundaryInput().get<std::string>(
                "Boundaries." + Field<T>::name() + "." + patch.name() + ".type",
                "");

        if (typeStr.empty())
            continue;

        if (typeStr == "fixed")
            boundaryType = FIXED;
        else if (typeStr == "normal_gradient")
            boundaryType = NORMAL_GRADIENT;
        else if (typeStr == "symmetry")
            boundaryType = SYMMETRY;
        else if (typeStr == "outflow")
            boundaryType = OUTFLOW;
        else
            throw Exception("FiniteVolumeField<T>", "setBoundaryTypes", "invalid boundary type \"" + typeStr + "\".");

        patchBoundaries_[patch.id()] = std::make_pair(boundaryType, T());
    }
}

//- External operators

template<class T>
FiniteVolumeField<T> operator+(FiniteVolumeField<T> lhs, const FiniteVolumeField<T> &rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
FiniteVolumeField<T> operator-(FiniteVolumeField<T> lhs, const FiniteVolumeField<T> &rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class T>
FiniteVolumeField<T> operator*(FiniteVolumeField<T> lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

template<class T>
FiniteVolumeField<T> operator*(Scalar lhs, FiniteVolumeField<T> rhs)
{
    rhs *= lhs;
    return rhs;
}

template<class T>
FiniteVolumeField<T> operator/(FiniteVolumeField<T> lhs, Scalar rhs)
{
    lhs /= rhs;
    return lhs;
}

//- External functions

template<class T, class TFunc>
void smooth(const FiniteVolumeField<T> &field,
            const CellGroup &cellsToSmooth,
            const CellGroup &cells,
            Scalar epsilon,
            FiniteVolumeField<T> &smoothedField,
            const TFunc &kernel)
{
    for (const Cell &cell: cellsToSmooth)
    {
        //- Determine the normalizing constant for this kernel
        Scalar integralK = 0.;

        auto kCells = cells.itemsWithin(Circle(cell.centroid(), epsilon));

        for (const Cell &kCell: kCells)
            integralK += kernel(cell, kCell, epsilon) * kCell.volume();

        Scalar tilde = 0.;
        for (const Cell &kCell: kCells)
            tilde += field(kCell) * kernel(cell, kCell, epsilon) * kCell.volume();

        smoothedField(cell) = tilde / integralK;
    }
}
