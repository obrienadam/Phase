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
    }
}

Scalar ScalarFiniteVolumeField::max() const
{
    const auto &self = *this;
    Scalar max = self[0];

    for(int i = 1, size = grid.cells.size(); i < size; ++i)
        max = std::max(self[i], max);

    for(int i = 0, size = grid.faces.size(); i < size; ++i)
        max = std::max(self.faces_[i], max);

    return max;
}

Scalar ScalarFiniteVolumeField::min() const
{
    const auto &self = *this;
    Scalar min = self[0];

    for(int i = 1, size = grid.cells.size(); i < size; ++i)
        min = std::min(self[i], min);

    for(int i = 0, size = grid.faces.size(); i < size; ++i)
        min = std::min(self.faces_[i], min);

    return min;
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
