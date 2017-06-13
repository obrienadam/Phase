#include "FiniteVolumeField.h"
#include "Exception.h"

#include <boost/algorithm/string.hpp>

//- Constructors

template<class T>
FiniteVolumeField<T>::FiniteVolumeField(const std::shared_ptr<const FiniteVolumeGrid2D>& grid,
                                        const std::string &name,
                                        const T& val,
                                        bool faces,
                                        bool nodes)
        :
        Field<T>::Field(grid->cells().size(), val, name),
        grid_(grid)
{
    if (faces)
        faces_.resize(grid_->faces().size(), val);

    if (nodes)
        nodes_.resize(grid_->nodes().size(), val);
}

template<class T>
FiniteVolumeField<T>::FiniteVolumeField(const Input &input,
                                        const std::shared_ptr<const FiniteVolumeGrid2D>& grid,
                                        const std::string &name,
                                        const T& val,
                                        bool faces,
                                        bool nodes)
        :
        FiniteVolumeField(grid, name, val, faces, nodes)
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

    if(!faces_.empty())
    {
        for (const Face &face: grid().interiorFaces())
            faces_[face.id()] = val;
    }

    if(!nodes_.empty())
    {
        for(const Node& node: grid().interiorNodes())
            nodes_[node.id()] = val;
    }
}

template<class T>
void FiniteVolumeField<T>::copyBoundaryTypes(const FiniteVolumeField &other)
{
    patchBoundaries_ = other.patchBoundaries_;
}

template<class T>
typename FiniteVolumeField<T>::BoundaryType FiniteVolumeField<T>::boundaryType(const Patch &patch) const
{
    return patchBoundaries_.find(patch.id())->second.first;
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
void FiniteVolumeField<T>::interpolateFaces(const std::function<Scalar(const Face &)> &alpha)
{
    auto &self = *this;

    for(const Face& face: grid().interiorFaces())
    {
        Scalar g = alpha(face);
        self(face) = g*self(face.lCell()) + (1. - g)*self(face.rCell());
    }

    setBoundaryFaces();
}

template<class T>
void FiniteVolumeField<T>::setBoundaryFaces()
{
    auto &self = *this;

    for(const Patch& patch: grid().patches())
    {
        switch (boundaryType(patch))
        {
            case FIXED:
                break;
            case NORMAL_GRADIENT:
            case SYMMETRY:
                for(const Face& face: patch)
                    self(face) = self(face.lCell());
                break;
        }
    }
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::savePreviousTimeStep(Scalar timeStep, int nPreviousFields)
{
    auto prevTimeStep = std::make_shared<PreviousField>(timeStep, *this);
    prevTimeStep->second.clearHistory();

    previousTimeSteps_.insert(previousTimeSteps_.begin(), prevTimeStep);
    previousTimeSteps_.resize(nPreviousFields);

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

template<class T>
Vector FiniteVolumeField<T>::vectorize() const
{
    const auto &self = *this;
    Vector vec = Vector(grid().nLocalActiveCells(), 0.);

    for (const Cell &cell: grid().localActiveCells())
        vec[cell.index(0)] = self(cell);

    return vec;
}

//- Operators

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator+=(const FiniteVolumeField &rhs)
{
    auto &self = *this;

    for(const Cell& cell: grid().cells())
        self(cell) += rhs(cell);

    if(!faces_.empty() && rhs.hasFaces())
        for(const Face& face: grid().faces())
            self(face) += rhs(face);

    if(!nodes_.empty() && rhs.hasNodes())
        for(const Node& node: grid().nodes())
            self(node) += rhs(node);

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator-=(const FiniteVolumeField &rhs)
{
    auto &self = *this;

    for(const Cell& cell: grid().cells())
        self(cell) -= rhs(cell);

    if(!faces_.empty() && rhs.hasFaces())
        for(const Face& face: grid().faces())
            self(face) -= rhs(face);

    if(!nodes_.empty() && rhs.hasNodes())
        for(const Node& node: grid().nodes())
            self(node) -= rhs(node);

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator*=(const FiniteVolumeField<Scalar> &rhs)
{
    auto &self = *this;

    for(const Cell& cell: grid().cells())
        self(cell) *= rhs(cell);

    if(!faces_.empty() && rhs.hasFaces())
        for(const Face& face: grid().faces())
            self(face) *= rhs(face);

    if(!nodes_.empty() && rhs.hasNodes())
        for(const Node& node: grid().nodes())
            self(node) *= rhs(node);

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator/=(const FiniteVolumeField<Scalar> &rhs)
{
    auto &self = *this;

    for(const Cell& cell: grid().cells())
        self(cell) /= rhs(cell);

    if(!faces_.empty() && rhs.hasFaces())
        for(const Face& face: grid().faces())
            self(face) /= rhs(face);

    if(!nodes_.empty() && rhs.hasNodes())
        for(const Node& node: grid().nodes())
            self(node) /= rhs(node);

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator*=(Scalar rhs)
{
    auto &self = *this;

    for(const Cell& cell: grid().cells())
        self(cell) *= rhs;

    if(!faces_.empty())
        for(const Face& face: grid().faces())
            self(face) *= rhs;

    if(!nodes_.empty())
        for(const Node& node: grid().nodes())
            self(node) *= rhs;

    return self;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator/=(Scalar rhs)
{
    auto &self = *this;

    for(const Cell& cell: grid().cells())
        self(cell) /= rhs;

    if(!faces_.empty())
        for(const Face& face: grid().faces())
            self(face) /= rhs;

    if(!nodes_.empty())
        for(const Node& node: grid().nodes())
            self(node) /= rhs;

    return self;
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

        for (const Patch& patch: grid().patches())
        {
            patchBoundaries_[patch.id()] = std::make_pair(boundaryType, T());
        }
    }

    for (const Patch& patch: grid().patches())
    {
        typeStr = input.boundaryInput().get<std::string>("Boundaries." + Field<T>::name() + "." + patch.name() + ".type",
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

template<class T>
void interpolateNodes(FiniteVolumeField<T> &field)
{
    for (const Node &node: field.grid().nodes())
    {
        const Scalar nCells = node.cells().size();
        field(node) = T();

        for (const Cell &cell: node.cells())
            field(node) += field(cell) / nCells;
    }

    for (const Face &face: field.grid().interiorFaces())
        field(face) = 0.5 * (field(face.lNode()) + field(face.rNode()));

    for(const Patch& patch: field.grid().patches())
    {
        switch (field.boundaryType(patch))
        {
            case FiniteVolumeField<T>::FIXED:
                break;

            case FiniteVolumeField<T>::NORMAL_GRADIENT:
            case FiniteVolumeField<T>::SYMMETRY:
                for(const Face& face: patch)
                    field(face) = field(face.lCell());
                break;

            default:
                throw Exception("FiniteVolumeField<T>", "interpolateNodes", "unrecongnized boundary condition type.");
        }
    }
}

template<class T>
FiniteVolumeField<T>
smooth(const FiniteVolumeField<T> &field, const std::vector<std::vector<Ref<const Cell> > > &rangeSearch, Scalar e)
{
    FiniteVolumeField<T> smoothedField(field.gridPtr(), field.name());
    Scalar A = 1.;
    const Scalar eSqr = e * e;

//    auto K = [&A](Scalar rSqr, Scalar eSqr){ return rSqr < eSqr ? A*pow(eSqr - rSqr, 3) : 0.; };
    auto K = [&A](Scalar rSqr, Scalar eSqr) { // This smoothing kernel appears to be slightly better
        Scalar r = sqrt(rSqr), e = sqrt(eSqr);

        return r < e ? A / (2. * e) * (1. + cos(M_PI * r / e)) : 0.;
    };

    for (const Cell &cell: field.grid().localActiveCells())
    {
        //- Determine the normalizing constant for this kernel
        Scalar integralK = 0.;
        A = 1.;

        for (const Cell &kCell: rangeSearch[cell.id()])
            integralK += K((kCell.centroid() - cell.centroid()).magSqr(), eSqr) * kCell.volume();

        A = 1. / integralK;

        smoothedField[cell.id()] = 0.;
        for (const Cell &kCell: rangeSearch[cell.id()])
            smoothedField[cell.id()] +=
                    field[kCell.id()] * K((kCell.centroid() - cell.centroid()).magSqr(), eSqr) * kCell.volume();
    }

    return smoothedField;
}
