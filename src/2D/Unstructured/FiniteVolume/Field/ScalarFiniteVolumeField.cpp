#include "ScalarFiniteVolumeField.h"

//- Protected methods

template<>
void ScalarFiniteVolumeField::setBoundaryRefValues(const Input &input)
{
    using namespace std;

    if (input.boundaryInput().count("Boundaries." + name_ + ".*") != 0)
    {
        Scalar refVal = input.boundaryInput().get<Scalar>("Boundaries." + name_ + ".*.value");

        for (const FaceGroup &patch: grid()->patches())
        {
            BoundaryType type = patchBoundaries_[patch.name()].first;
            patchBoundaries_[patch.name()] = std::make_pair(type, refVal);
        }
    }

    for (const FaceGroup &patch: grid()->patches())
    {
        Scalar refVal = input.boundaryInput().get<Scalar>("Boundaries." + name_ + "." + patch.name() + ".value", 0);
        BoundaryType type = patchBoundaries_[patch.name()].first;
        patchBoundaries_[patch.name()] = std::make_pair(type, refVal);
    }

    auto &self = *this;

    for (const FaceGroup &patch: grid()->patches())
        for (const Face &face: patch)
            self(face) = boundaryRefValue(patch);
}

template<>
bool ScalarFiniteVolumeField::isfinite() const
{
    int isFinite = 1;
    for(const Cell &cell: cells())
        if(!std::isfinite((*this)(cell)))
        {
            isFinite = 0;
            break;
        }

    return bool(grid_->comm().min(isFinite));
}

//- External functions

ScalarFiniteVolumeField operator*(const ScalarFiniteVolumeField &lhs, ScalarFiniteVolumeField rhs)
{
    rhs *= lhs;
    return rhs;
}

ScalarFiniteVolumeField operator-(const ScalarFiniteVolumeField &lhs, Scalar rhs)
{
    ScalarFiniteVolumeField diff(lhs.grid(), "", 0., false);

    std::transform(lhs.begin(), lhs.end(), diff.begin(), [rhs](Scalar val) {
        return val - rhs;
    });

    return lhs;
}

ScalarFiniteVolumeField operator/(ScalarFiniteVolumeField lhs, const ScalarFiniteVolumeField &rhs)
{
    lhs /= rhs;
    return lhs;
}

ScalarFiniteVolumeField operator/(Scalar lhs, ScalarFiniteVolumeField rhs)
{
    for (const Cell &cell: rhs.grid()->cells())
        rhs(cell) = lhs / rhs(cell);

    if (rhs.hasFaces())
        for (const Face &face: rhs.grid()->faces())
            rhs(face) = lhs / rhs(face);

    if (rhs.hasNodes())
        for (const Node &node: rhs.grid()->nodes())
            rhs(node) = lhs / rhs(node);

    return rhs;
}
