#include "FiniteVolumeField.h"
#include "Exception.h"

//- Constructors

template <class T>
FiniteVolumeField<T>::FiniteVolumeField(const FiniteVolumeGrid2D &grid, const std::string &name)
    :
      Field<T>::Field(grid.cells().size(), 0., name),
      faces_(grid.faces().size(), 0.),
      grid(grid),
      prevFieldPtr_(nullptr)
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
      boundaryTypes_(other.boundaryTypes_),
      boundaryRefValues_(other.boundaryRefValues_),
      faces_(other.faces_),
      prevFieldPtr_(nullptr)
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
    for(const Face &face: grid.interiorFaces())
        faces_[face.id()] = val;
}

template<class T>
void FiniteVolumeField<T>::copyBoundaryTypes(const FiniteVolumeField &other)
{
    boundaryTypes_ = other.boundaryTypes_;
    boundaryRefValues_.resize(other.boundaryRefValues_.size());
}

template<class T>
typename FiniteVolumeField<T>::BoundaryType FiniteVolumeField<T>::boundaryType(size_t faceId) const
{
    if(boundaryTypes_.size() == 0)
        return NORMAL_GRADIENT;

    return boundaryTypes_[grid.faces()[faceId].patch().id()];
}

template<class T>
T FiniteVolumeField<T>::boundaryRefValue(size_t faceId) const
{
    if(boundaryRefValues_.size() == 0)
        return T();

    return boundaryRefValues_[grid.faces()[faceId].patch().id()];
}

template<class T>
FiniteVolumeField<T>& FiniteVolumeField<T>::save()
{
    prevFieldPtr_ = std::shared_ptr<FiniteVolumeField>(new FiniteVolumeField<T>(*this));
    return *prevFieldPtr_;
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
    boundaryTypes_ = rhs.boundaryTypes_;
    boundaryRefValues_ = rhs.boundaryRefValues_;
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

//- Protected methods

template<class T>
void FiniteVolumeField<T>::setBoundaryTypes(const Input &input)
{
    using namespace std;

    //- Check for a default patch

    auto defBoundary = input.boundaryInput().get_child_optional("Boundaries." + Field<T>::name + ".*");

    if(defBoundary)
    {
        std::string type = input.boundaryInput().get<string>("Boundaries." + Field<T>::name + ".*.type");

        for(int i = 0; i < grid.patches().size(); ++i)
        {
            if(type == "fixed")
                boundaryTypes_.push_back(FIXED);
            else if(type == "normal_gradient")
                boundaryTypes_.push_back(NORMAL_GRADIENT);
            else
                throw Exception("FiniteVolumeField", "FiniteVolumeField", "unrecognized boundary type \"" + type + "\".");
        }
    }
    else
    {
        //- Boundary condition association to patches is done by patch id
        for(const Patch& patch: grid.patches())
        {
            string root = "Boundaries." + Field<T>::name + "." + patch.name;
            string type = input.boundaryInput().get<string>(root + ".type");

            if(type == "fixed")
                boundaryTypes_.push_back(FIXED);
            else if (type == "normal_gradient")
                boundaryTypes_.push_back(NORMAL_GRADIENT);
            else
                throw Exception("FiniteVolumeField", "FiniteVolumeField", "unrecognized boundary type \"" + type + "\".");
        }
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
        }
    }
}

template<class T>
FiniteVolumeField<T> smooth(const FiniteVolumeField<T>& field, const RangeSearch& rangeSearch, Scalar h)
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

        for(const Cell &kCell: rangeSearch.getResult(cell.id()))
        {
            totalVol += kCell.volume();
            intKr += kr((cell.centroid() - kCell.centroid()).mag(), h);
        }

        for(const Cell &kCell: rangeSearch.getResult(cell.id()))
        {
            smoothedField[cell.id()] += field[kCell.id()]*kr((cell.centroid() - kCell.centroid()).mag(), h)/intKr;
        }

        //smoothedField[cell.id()] = std::min(field[cell.id()], std::max(field[cell.id()], smoothedField[cell.id()]));
    }

    return smoothedField;
}
