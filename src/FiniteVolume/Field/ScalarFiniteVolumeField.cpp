#include "ScalarFiniteVolumeField.h"
#include "FiniteVolumeGrid2D.h"
#include "Exception.h"

ScalarFiniteVolumeField::ScalarFiniteVolumeField(const FiniteVolumeGrid2D &grid, const std::string &name)
    :
      Field<Scalar>::Field(grid.cells.size(), 0., name),
      faces_(grid.faces.size(), 0.),
      grid(grid)
{

}

ScalarFiniteVolumeField::ScalarFiniteVolumeField(const Input &input, const FiniteVolumeGrid2D &grid, const std::string &name)
    :
      ScalarFiniteVolumeField(grid, name)
{
    auto &self = *this;

    //- Boundary condition association to patches is done by patch id
    for(const Patch& patch: grid.patches())
    {
        std::string root = "Boundaries." + name + "." + patch.name;
        std::string type = input.boundaryInput().get<std::string>(root + ".type");

        if(type == "fixed")
        {
            boundaryTypes_.push_back(FIXED);
        }
        else if (type == "normal_gradient")
        {
            boundaryTypes_.push_back(NORMAL_GRADIENT);
        }
        else
        {
            throw Exception("ScalarFiniteVolumeField", "ScalarFiniteVolumeField", "unrecognized boundary type \"" + type + "\".");
        }

        boundaryRefValues_.push_back(input.boundaryInput().get<Scalar>(root + ".value"));

        for(const Face& face: patch.faces())
            self.faces()[face.id()] = boundaryRefValues_.back();
    }
}

void ScalarFiniteVolumeField::fill(Scalar val)
{
    std::fill(begin(), end(), val);
    std::fill(faces_.begin(), faces_.end(), val);
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
    return boundaryTypes_[grid.faces[faceId].patch().id()];
}

ScalarFiniteVolumeField& ScalarFiniteVolumeField::operator *=(const ScalarFiniteVolumeField& rhs)
{
    auto &self = *this;

    for(int i = 0, end = self.size(); i < end; ++i)
        self[i] *= rhs[i];

    return self;
}

//- External functions

ScalarFiniteVolumeField operator*(const ScalarFiniteVolumeField& lhs, ScalarFiniteVolumeField rhs)
{
    rhs *= lhs;
    return rhs;
}
