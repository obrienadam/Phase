#include "ScalarFiniteVolumeField.h"

template<>
ScalarFiniteVolumeField& ScalarFiniteVolumeField::operator =(const SparseVector& rhs)
{
    auto &self = *this;

    for(const Cell &cell: self.grid.activeCells())
        self[cell.id()] = rhs[cell.globalIndex()];

    return self;
}

//- Protected methods

template<>
void ScalarFiniteVolumeField::setBoundaryRefValues(const Input &input)
{
    using namespace std;

    if(input.boundaryInput().count("Boundaries." + name + ".*") != 0)
    {
        Scalar refVal = input.boundaryInput().get<Scalar>("Boundaries." + name + ".*.value");

        for(const auto &entry: grid.patches())
        {
            BoundaryType type = patchBoundaries_[entry.second.id()].first;
            patchBoundaries_[entry.second.id()] = std::make_pair(type, refVal);
        }
    }

    for(const auto &entry: grid.patches())
    {
        Scalar refVal = input.boundaryInput().get<Scalar>("Boundaries." + name + "." + entry.first + ".value", 0);
        BoundaryType type = patchBoundaries_[entry.second.id()].first;
        patchBoundaries_[entry.second.id()] = std::make_pair(type, refVal);
    }

    auto &self = *this;

    for(const Face &face: grid.boundaryFaces())
        self(face) = boundaryRefValue(face);
}

//- External functions

ScalarFiniteVolumeField operator*(const ScalarFiniteVolumeField& lhs, ScalarFiniteVolumeField rhs)
{
    rhs *= lhs;
    return rhs;
}

ScalarFiniteVolumeField operator/(ScalarFiniteVolumeField lhs, const ScalarFiniteVolumeField& rhs)
{
    lhs /= rhs;
    return lhs;
}
