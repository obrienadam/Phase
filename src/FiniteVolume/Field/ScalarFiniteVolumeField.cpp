#include "ScalarFiniteVolumeField.h"

template<>
ScalarFiniteVolumeField& ScalarFiniteVolumeField::operator =(const Vector& rhs)
{
    auto &self = *this;

    for(const Cell &cell: self.grid.localActiveCells())
        self[cell.id()] = rhs[cell.localIndex()];

    return self;
}

//- Protected methods

template<>
void ScalarFiniteVolumeField::setBoundaryRefValues(const Input &input)
{
    using namespace std;

    if(input.boundaryInput().count("Boundaries." + name_ + ".*") != 0)
    {
        Scalar refVal = input.boundaryInput().get<Scalar>("Boundaries." + name_ + ".*.value");

        for(const auto &entry: grid.patches())
        {
            BoundaryType type = patchBoundaries_[entry.second.id()].first;
            patchBoundaries_[entry.second.id()] = std::make_pair(type, refVal);
        }
    }

    for(const auto &entry: grid.patches())
    {
        Scalar refVal = input.boundaryInput().get<Scalar>("Boundaries." + name_ + "." + entry.first + ".value", 0);
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

ScalarFiniteVolumeField operator/(Scalar lhs, ScalarFiniteVolumeField rhs)
{
    for(const Cell& cell: rhs.grid.cells())
        rhs(cell) = lhs/rhs(cell);

    for(const Face& face: rhs.grid.faces())
        rhs(face) = lhs/rhs(face);

    return rhs;
}
