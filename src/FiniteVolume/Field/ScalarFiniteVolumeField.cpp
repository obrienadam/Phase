#include "ScalarFiniteVolumeField.h"

template<>
ScalarFiniteVolumeField &ScalarFiniteVolumeField::operator=(const Vector &rhs)
{
    auto &self = *this;

    for (const Cell &cell: self.grid().localActiveCells())
        self[cell.id()] = rhs[cell.index(0)];

    return self;
}

//- Protected methods

template<>
void ScalarFiniteVolumeField::setBoundaryRefValues(const Input &input)
{
    using namespace std;

    if (input.boundaryInput().count("Boundaries." + name_ + ".*") != 0)
    {
        Scalar refVal = input.boundaryInput().get<Scalar>("Boundaries." + name_ + ".*.value");

        for (const Patch& patch: grid().patches())
        {
            BoundaryType type = patchBoundaries_[patch.id()].first;
            patchBoundaries_[patch.id()] = std::make_pair(type, refVal);
        }
    }

    for (const Patch& patch: grid().patches())
    {
        Scalar refVal = input.boundaryInput().get<Scalar>("Boundaries." + name_ + "." + patch.name() + ".value", 0);
        BoundaryType type = patchBoundaries_[patch.id()].first;
        patchBoundaries_[patch.id()] = std::make_pair(type, refVal);
    }

    auto &self = *this;

    for(const Patch& patch: grid().patches())
        for(const Face& face: patch)
            self(face) = boundaryRefValue(patch);
}

//- External functions

ScalarFiniteVolumeField operator*(const ScalarFiniteVolumeField &lhs, ScalarFiniteVolumeField rhs)
{
    rhs *= lhs;
    return rhs;
}

ScalarFiniteVolumeField operator/(ScalarFiniteVolumeField lhs, const ScalarFiniteVolumeField &rhs)
{
    lhs /= rhs;
    return lhs;
}

ScalarFiniteVolumeField operator/(Scalar lhs, ScalarFiniteVolumeField rhs)
{
    for (const Cell &cell: rhs.grid().cells())
        rhs(cell) = lhs / rhs(cell);

    for (const Face &face: rhs.grid().faces())
        rhs(face) = lhs / rhs(face);

    return rhs;
}
