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
      previousFields_()
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
    interpolateFaces(*this);
}

template<class T>
void FiniteVolumeField<T>::copyBoundaryTypes(const FiniteVolumeField &other)
{
    patchBoundaries_ = other.patchBoundaries_;
}

template<class T>
typename FiniteVolumeField<T>::BoundaryType FiniteVolumeField<T>::boundaryType(size_t faceId) const
{
    if(patchBoundaries_.size() == 0)
        return NORMAL_GRADIENT;

    return (patchBoundaries_.find(grid.faces()[faceId].patch().id())->second).first;
}

template<class T>
T FiniteVolumeField<T>::boundaryRefValue(size_t faceId) const
{
    if(patchBoundaries_.size() == 0)
        return T();

    return (patchBoundaries_.find(grid.faces()[faceId].patch().id())->second).second;
}

template<class T>
std::pair<typename FiniteVolumeField<T>::BoundaryType, T> FiniteVolumeField<T>::boundaryInfo(size_t faceId) const
{
    if(patchBoundaries_.size() == 0)
        return std::make_pair(NORMAL_GRADIENT, T());

    return patchBoundaries_.find(grid.faces()[faceId].patch().id())->second;
}

template<class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::save(Scalar timeStep, int nPreviousFields)
{
    previousFields_.push_front(std::make_pair(timeStep, FiniteVolumeField<T>(*this)));

    while(previousFields_.size() > nPreviousFields)
        previousFields_.pop_back();

    return previousFields_.front().second;
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


    std::string typeStr = input.boundaryInput().get<string>("Boundaries." + Field<T>::name + ".*.type", "");
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
        else
            throw Exception("FiniteVolumeField<T>", "setBoundaryTypes", "invalid boundary type \"" + typeStr + "\".");

        for(const auto &entry: grid.patches())
        {
            patchBoundaries_[entry.second.id()] = std::make_pair(boundaryType, T());
        }
    }

    for(const auto &entry: grid.patches())
    {
        typeStr = input.boundaryInput().get<std::string>("Boundaries." + Field<T>::name + "." + entry.first + ".type", "");

        if(typeStr.empty())
            continue;

        if(typeStr == "fixed")
            boundaryType = FIXED;
        else if(typeStr == "normal_gradient")
            boundaryType = NORMAL_GRADIENT;
        else if(typeStr == "symmetry")
            boundaryType = SYMMETRY;
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
void interpolateFaces(FiniteVolumeField<T>& field)
{
    Vector2D rf, sf;

    for(const Face& face: field.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        Scalar alpha = rCell.volume()/(lCell.volume() + rCell.volume());
        field.faces()[face.id()] = field[lCell.id()]*alpha + field[rCell.id()]*(1. - alpha);
    }

    for(const Face& face: field.grid.boundaryFaces())
    {
        switch(field.boundaryType(face.id()))
        {
        case FiniteVolumeField<T>::FIXED:
            break;

        case FiniteVolumeField<T>::NORMAL_GRADIENT:
            rf = face.centroid() - face.lCell().centroid();
            sf = face.outwardNorm(face.lCell().centroid());
            field.faces()[face.id()] = sf.mag()/(dot(rf, sf)/dot(rf, rf))*field.boundaryRefValue(face.id()) + field[face.lCell().id()];
            break;

        case FiniteVolumeField<T>::SYMMETRY:
            field.faces()[face.id()] = field[face.lCell().id()];
            break;

        default:
            throw Exception("FiniteVolumeField<T>", "interpolateFaces", "unrecongnized boundary condition type.");
        }
    }
}

template<class T>
void harmonicInterpolateFaces(FiniteVolumeField<T>& field)
{
    Vector2D rf, sf;

    for(const Face& face: field.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        Scalar alpha = rCell.volume()/(lCell.volume() + rCell.volume());
        field.faces()[face.id()] = 1./(alpha/field[lCell.id()] + (1. - alpha)/field[rCell.id()]);
    }

    for(const Face& face: field.grid.boundaryFaces())
    {
        switch(field.boundaryType(face.id()))
        {
        case FiniteVolumeField<T>::FIXED:
            break;

        case FiniteVolumeField<T>::NORMAL_GRADIENT:
            rf = face.centroid() - face.lCell().centroid();
            sf = face.outwardNorm(face.lCell().centroid());
            field.faces()[face.id()] = sf.mag()/(dot(rf, sf)/dot(rf, rf))*field.boundaryRefValue(face.id()) + field[face.lCell().id()];
            break;

        case FiniteVolumeField<T>::SYMMETRY:
            field.faces()[face.id()] = field[face.lCell().id()];
            break;

        default:
            throw Exception("FiniteVolumeField<T>", "harmonicInterpolateFaces", "unrecongnized boundary condition type.");
        }
    }
}

template<class T>
FiniteVolumeField<T> smooth(const FiniteVolumeField<T>& field, const std::vector< std::vector< Ref<const Cell> > >& rangeSearch, Scalar h)
{
    FiniteVolumeField<T> smoothedField(field.grid, field.name);

    auto pow4 = [](Scalar x) { return x*x*x*x; };
    auto pow2 = [](Scalar x) { return x*x; };
    auto kr = [&pow2, &pow4](Scalar r, Scalar h)
    {
        return pow4(1. - pow2(std::min(r/h, 1.)));
    };

    for(const Cell &cell: field.grid.cells())
    {
        Scalar totalVol = 0., intKr = 0.;

        for(const Cell &kCell: rangeSearch[cell.id()])
        {
            totalVol += kCell.volume();
            intKr += kr((cell.centroid() - kCell.centroid()).mag(), h);
        }

        for(const Cell &kCell: rangeSearch[cell.id()])
        {
            smoothedField[cell.id()] += field[kCell.id()]*kr((cell.centroid() - kCell.centroid()).mag(), h)/intKr;
        }

        //smoothedField[cell.id()] = std::min(field[cell.id()], std::max(field[cell.id()], smoothedField[cell.id()]));
    }

    return smoothedField;
}
