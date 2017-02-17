#include "FiniteVolumeField.h"
#include "Exception.h"

#include <boost/algorithm/string.hpp>

//- Constructors

template <class T>
FiniteVolumeField<T>::FiniteVolumeField(const FiniteVolumeGrid2D &grid, const std::string &name)
    :
      Field<T>::Field(grid.cells().size(), 0., name),
      faces_(grid.faces().size(), 0.),
      grid(grid)
{

}

template<class T>
FiniteVolumeField<T>::FiniteVolumeField(const Input &input, const FiniteVolumeGrid2D &grid, const std::string &name)
    :
      FiniteVolumeField(grid, name)
{
    using namespace std;
    setBoundaryTypes(input);
    setBoundaryRefValues(input);
}

template<class T>
FiniteVolumeField<T>::FiniteVolumeField(const FiniteVolumeField &other)
    :
      Field<T>::Field(other),
      grid(other.grid),
      patchBoundaries_(other.patchBoundaries_),
      faces_(other.faces_),
      previousTimeSteps_(),
      previousIteration_()
{

}

//- Public methods

template<class T>
void FiniteVolumeField<T>::fill(const T &val)
{
    std::fill(FiniteVolumeField<T>::begin(), FiniteVolumeField<T>::end(), val);
    std::fill(faces_.begin(), faces_.end(), val);
}

template<class T>
void FiniteVolumeField<T>::fillInterior(const T &val)
{
    std::fill(FiniteVolumeField<T>::begin(), FiniteVolumeField<T>::end(), val);
    for(const Face& face: grid.interiorFaces())
        faces_[face.id()] = val;
}

template<class T>
void FiniteVolumeField<T>::copyBoundaryTypes(const FiniteVolumeField &other)
{
    patchBoundaries_ = other.patchBoundaries_;
}

template<class T>
typename FiniteVolumeField<T>::BoundaryType FiniteVolumeField<T>::boundaryType(const Face &face) const
{
    if(patchBoundaries_.size() == 0)
        return NORMAL_GRADIENT;

    return (patchBoundaries_.find(grid.faces()[face.id()].patch().id())->second).first;
}

template<class T>
T FiniteVolumeField<T>::boundaryRefValue(const Face &face) const
{
    if(patchBoundaries_.size() == 0)
        return T();

    return (patchBoundaries_.find(grid.faces()[face.id()].patch().id())->second).second;
}

template<class T>
std::pair<typename FiniteVolumeField<T>::BoundaryType, T> FiniteVolumeField<T>::boundaryInfo(const Face &face) const
{
    if(patchBoundaries_.size() == 0)
        return std::make_pair(NORMAL_GRADIENT, T());

    return patchBoundaries_.find(grid.faces()[face.id()].patch().id())->second;
}

template<class T>
void FiniteVolumeField<T>::setBoundaryFaces()
{
    for(const Face &face: grid.boundaryFaces())
    {
        switch(boundaryType(face))
        {
        case FIXED:
            break;
        case NORMAL_GRADIENT: case SYMMETRY:
            faces_[face.id()] = std::vector<T>::operator [](face.lCell().id());
            break;
        }
    }
}

template<class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::savePreviousTimeStep(Scalar timeStep, int nPreviousFields)
{
    previousTimeSteps_.push_front(std::make_pair(timeStep, FiniteVolumeField<T>(*this)));

    while(previousTimeSteps_.size() > nPreviousFields)
        previousTimeSteps_.pop_back();

    return previousTimeSteps_.front().second;
}

template<class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::savePreviousIteration()
{
    if(previousIteration_.size() >= 1)
        previousIteration_.clear();

    previousIteration_.push_back(FiniteVolumeField<T>(*this));

    return previousIteration_.front();
}

template<class T>
Vector FiniteVolumeField<T>::vectorize() const
{
    const auto& self = *this;
    Vector vec = Vector(grid.nLocalActiveCells(), 0.);

    for(const Cell& cell: grid.localActiveCells())
        vec[cell.index(0)] = self(cell);

    return vec;
}

//- Operators

template<class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::operator=(const FiniteVolumeField& rhs)
{
    if(this == &rhs)
        return *this;
    else if(&grid != &rhs.grid)
        throw Exception("FiniteVolumeField", "operator=", "grid references must be the same.");

    Field<T>::operator =(rhs);
    patchBoundaries_ = rhs.patchBoundaries_;
    faces_ = rhs.faces_;

    return *this;
}

template <class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::operator+=(const FiniteVolumeField& rhs)
{
    auto &self = *this;

    for(size_t i = 0, end = self.size(); i < end; ++i)
        self[i] += rhs[i];

    for(size_t i = 0, end = self.faces().size(); i < end; ++i)
        self.faces()[i] += rhs.faces()[i];

    return *this;
}

template <class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::operator-=(const FiniteVolumeField& rhs)
{
    auto &self = *this;

    for(size_t i = 0, end = self.size(); i < end; ++i)
        self[i] -= rhs[i];

    for(size_t i = 0, end = self.faces().size(); i < end; ++i)
        self.faces()[i] -= rhs.faces()[i];

    return *this;
}

template<class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::operator*=(const FiniteVolumeField<Scalar>& rhs)
{
    auto &self = *this;

    for(int i = 0, end = self.size(); i < end; ++i)
        self[i] *= rhs[i];

    for(int i = 0, end = self.faces().size(); i < end; ++i)
        self.faces()[i] *= rhs.faces()[i];

    return self;
}

template <class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::operator*=(Scalar rhs)
{
    auto &self = *this;

    for(size_t i = 0, end = self.size(); i < end; ++i)
        self[i] *= rhs;

    for(size_t i = 0, end = self.faces().size(); i < end; ++i)
        self.faces()[i] *= rhs;

    return *this;
}

template <class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::operator/=(Scalar rhs)
{
    auto &self = *this;

    for(size_t i = 0, end = self.size(); i < end; ++i)
        self[i] /= rhs;

    for(size_t i = 0, end = self.faces().size(); i < end; ++i)
        self.faces()[i] /= rhs;

    return *this;
}

template<class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::operator/=(const FiniteVolumeField<Scalar>& rhs)
{
    auto &self = *this;

    for(int i = 0, end = self.size(); i < end; ++i)
        self[i] /= rhs[i];

    for(int i = 0, end = self.faces().size(); i < end; ++i)
        self.faces()[i] /= rhs.faces()[i];

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

    if(!typeStr.empty())
    {
        if(typeStr == "fixed")
            boundaryType = FIXED;
        else if(typeStr == "normal_gradient")
            boundaryType = NORMAL_GRADIENT;
        else if(typeStr == "symmetry")
            boundaryType = SYMMETRY;
        else if(typeStr == "outflow")
            boundaryType = OUTFLOW;
        else
            throw Exception("FiniteVolumeField<T>", "setBoundaryTypes", "invalid boundary type \"" + typeStr + "\".");

        for(const auto &entry: grid.patches())
        {
            patchBoundaries_[entry.second.id()] = std::make_pair(boundaryType, T());
        }
    }

    for(const auto &entry: grid.patches())
    {
        typeStr = input.boundaryInput().get<std::string>("Boundaries." + Field<T>::name() + "." + entry.first + ".type", "");

        if(typeStr.empty())
            continue;

        if(typeStr == "fixed")
            boundaryType = FIXED;
        else if(typeStr == "normal_gradient")
            boundaryType = NORMAL_GRADIENT;
        else if(typeStr == "symmetry")
            boundaryType = SYMMETRY;
        else if(typeStr == "outflow")
            boundaryType = OUTFLOW;
        else
            throw Exception("FiniteVolumeField<T>", "setBoundaryTypes", "invalid boundary type \"" + typeStr + "\".");

        patchBoundaries_[entry.second.id()] = std::make_pair(boundaryType, T());
    }
}

//- External operators

template<class T>
FiniteVolumeField<T> operator+(FiniteVolumeField<T> lhs, const FiniteVolumeField<T>& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
FiniteVolumeField<T> operator-(FiniteVolumeField<T> lhs, const FiniteVolumeField<T>& rhs)
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
    for(const Node& node: field.grid.nodes())
    {
        const Scalar nCells = node.cells().size();
        field(node) = T();

        for(const Cell &cell: node.cells())
            field(node) += field(cell)/nCells;
    }

    for(const Face& face: field.grid.interiorFaces())
        field(face) = 0.5*(field(face.lNode()) + field(face.rNode()));

    for(const Face& face: field.grid.boundaryFaces())
    {
        switch(field.boundaryType(face))
        {
        case FiniteVolumeField<T>::FIXED:
            break;

        case FiniteVolumeField<T>::NORMAL_GRADIENT: case FiniteVolumeField<T>::SYMMETRY:
            field(face) = field(face.lCell());
            break;

        default:
            throw Exception("FiniteVolumeField<T>", "interpolateNodes", "unrecongnized boundary condition type.");
        }
    }
}

template<class T>
FiniteVolumeField<T> smooth(const FiniteVolumeField<T>& field, const std::vector< std::vector< Ref<const Cell> > >& rangeSearch, Scalar e)
{
    FiniteVolumeField<T> smoothedField(field.grid, field.name());
    Scalar A = 1.;
    const Scalar eSqr = e*e;

//    auto K = [&A](Scalar rSqr, Scalar eSqr){ return rSqr < eSqr ? A*pow(eSqr - rSqr, 3) : 0.; };
    auto K = [&A](Scalar rSqr, Scalar eSqr){ // This smoothing kernel appears to be slightly better
        Scalar r = sqrt(rSqr), e = sqrt(eSqr);

        return r < e ? A/(2.*e)*(1. + cos(M_PI*r/e)) : 0.;
    };

    for(const Cell &cell: field.grid.localActiveCells())
    {
        //- Determine the normalizing constant for this kernel
        Scalar integralK = 0.;
        A = 1.;

        for(const Cell &kCell: rangeSearch[cell.id()])
            integralK += K((kCell.centroid() - cell.centroid()).magSqr(), eSqr)*kCell.volume();

        A = 1./integralK;

        smoothedField[cell.id()] = 0.;
        for(const Cell &kCell: rangeSearch[cell.id()])
            smoothedField[cell.id()] += field[kCell.id()]*K((kCell.centroid() - cell.centroid()).magSqr(), eSqr)*kCell.volume();
    }

    return smoothedField;
}
