#include "ScalarFiniteVolumeField.h"
#include "FiniteVolumeGrid2D.h"
#include "Exception.h"

ScalarFiniteVolumeField::ScalarFiniteVolumeField(const FiniteVolumeGrid2D &grid, const std::string &name)
    :
      Field<Scalar>::Field(grid.cells().size(), 0., name),
      faces_(grid.faces().size(), 0.),
      grid(grid)
{

}

ScalarFiniteVolumeField::ScalarFiniteVolumeField(const Input &input, const FiniteVolumeGrid2D &grid, const std::string &name)
    :
      ScalarFiniteVolumeField(grid, name)
{
    using namespace std;

    auto &self = *this;

    //- Check for a default patch

    auto defBoundary = input.boundaryInput().get_child_optional("Boundaries." + name + ".*");

    if(defBoundary)
    {
        std::string type = input.boundaryInput().get<string>("Boundaries." + name + ".*.type");

        for(int i = 0; i < grid.patches().size(); ++i)
        {
            if(type == "fixed")
                boundaryTypes_.push_back(FIXED);
            else if(type == "normal_gradient")
                boundaryTypes_.push_back(NORMAL_GRADIENT);
            else
                throw Exception("ScalarFiniteVolumeField", "ScalarFiniteVolumeField", "unrecognized boundary type \"" + type + "\".");

            boundaryRefValues_.push_back(input.boundaryInput().get<Scalar>("Boundaries." + name + ".*.value"));
        }
    }
    else
    {
        //- Boundary condition association to patches is done by patch id
        for(const Patch& patch: grid.patches())
        {
            string root = "Boundaries." + name + "." + patch.name;
            string type = input.boundaryInput().get<string>(root + ".type");

            if(type == "fixed")
                boundaryTypes_.push_back(FIXED);
            else if (type == "normal_gradient")
                boundaryTypes_.push_back(NORMAL_GRADIENT);
            else
                throw Exception("ScalarFiniteVolumeField", "ScalarFiniteVolumeField", "unrecognized boundary type \"" + type + "\".");

            boundaryRefValues_.push_back(input.boundaryInput().get<Scalar>(root + ".value"));
        }
    }

    for(const Patch& patch: grid.patches())
    {
        for(const Face& face: patch.faces())
            self.faces()[face.id()] = boundaryRefValue(face.id());
    }
}

void ScalarFiniteVolumeField::fill(Scalar val)
{
    std::fill(begin(), end(), val);
    std::fill(faces_.begin(), faces_.end(), val);
}

void ScalarFiniteVolumeField::fillInterior(Scalar val)
{
    std::fill(begin(), end(), val);

    for(const Face &face: grid.interiorFaces())
        faces_[face.id()] = val;
}

void ScalarFiniteVolumeField::copyBoundaryTypes(const ScalarFiniteVolumeField &other)
{
    boundaryTypes_ = other.boundaryTypes_;
    boundaryRefValues_.resize(other.boundaryRefValues_.size());
}

ScalarFiniteVolumeField& ScalarFiniteVolumeField::operator =(const SparseVector& rhs)
{
    auto &self = *this;

    for(int i = 0, end = self.size(); i < end; ++i)
        self[i] = rhs[i];

    return self;
}

ScalarFiniteVolumeField::BoundaryType ScalarFiniteVolumeField::boundaryType(size_t faceId) const
{
    if(boundaryTypes_.size() == 0)
        return ScalarFiniteVolumeField::NORMAL_GRADIENT;

    return boundaryTypes_[grid.faces()[faceId].patch().id()];
}

Scalar ScalarFiniteVolumeField::boundaryRefValue(size_t faceId) const
{
    if(boundaryRefValues_.size() == 0)
        return 0;

    return boundaryRefValues_[grid.faces()[faceId].patch().id()];
}

ScalarFiniteVolumeField& ScalarFiniteVolumeField::operator *=(const ScalarFiniteVolumeField& rhs)
{
    auto &self = *this;

    for(int i = 0, end = self.size(); i < end; ++i)
        self[i] *= rhs[i];

    for(int i = 0, end = self.faces().size(); i < end; ++i)
        self.faces()[i] *= rhs.faces()[i];

    return self;
}

//- External functions

ScalarFiniteVolumeField operator*(const ScalarFiniteVolumeField& lhs, ScalarFiniteVolumeField rhs)
{
    rhs *= lhs;
    return rhs;
}

void interpolateFaces(ScalarFiniteVolumeField& field)
{
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
        case ScalarFiniteVolumeField::FIXED:
            break;

        case ScalarFiniteVolumeField::NORMAL_GRADIENT:
            field.faces()[face.id()] = field[face.lCell().id()];
            break;
        }
    }
}
