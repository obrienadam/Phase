#include "VectorFiniteVolumeField.h"
#include "Exception.h"

template<>
VectorFiniteVolumeField& VectorFiniteVolumeField::operator=(const Vector& rhs)
{
    auto &self = *this;
    const size_t nActiveCells = grid.nLocalActiveCells();

    for(const Cell &cell: grid.localActiveCells())
    {
        self[cell.id()].x = rhs[cell.index(0)];
        self[cell.id()].y = rhs[cell.index(0) + nActiveCells];
    }

    return self;
}

template<>
Vector VectorFiniteVolumeField::vectorize() const
{
    const auto& self = *this;
    const Size nActiveCells = grid.nLocalActiveCells();
    Vector vec(2*nActiveCells, 0.);

    for(const Cell& cell: grid.localActiveCells())
    {
        vec[cell.index(0)] = self(cell).x;
        vec[nActiveCells + cell.index(0)] = self(cell).y;
    }

    return vec;
}

//- Protected methods

template<>
void VectorFiniteVolumeField::setBoundaryRefValues(const Input &input)
{
    using namespace std;

    string valStr = input.boundaryInput().get<string>("Boundaries." + name_ + ".*.value", "");

    if(!valStr.empty())
    {
        Vector2D refVal = Vector2D(valStr);

        for(const auto &entry: grid.patches())
        {
            BoundaryType type = patchBoundaries_[entry.second.id()].first;
            patchBoundaries_[entry.second.id()] = std::make_pair(type, refVal);
        }
    }

    for(const auto &entry: grid.patches())
    {
        valStr = input.boundaryInput().get<string>("Boundaries." + name_ + "." + entry.first + ".value", "");

        if(valStr.empty())
            continue;

        Vector2D refVal(valStr);

        BoundaryType type = patchBoundaries_[entry.second.id()].first;
        patchBoundaries_[entry.second.id()] = std::make_pair(type, refVal);
    }

    auto &self = *this;

    for(const Face &face: grid.boundaryFaces())
        self(face) = boundaryRefValue(face);
}

//- External functions

VectorFiniteVolumeField operator*(const ScalarFiniteVolumeField& lhs, VectorFiniteVolumeField rhs)
{
    rhs *= lhs;
    return rhs;
}

VectorFiniteVolumeField operator*(VectorFiniteVolumeField lhs, const ScalarFiniteVolumeField& rhs)
{
    lhs *= rhs;
    return lhs;
}

VectorFiniteVolumeField operator*(const ScalarFiniteVolumeField& lhs, const Vector2D& rhs)
{
    VectorFiniteVolumeField result(lhs.grid, lhs.name());

    for(const Cell& cell: lhs.grid.cells())
        result[cell.id()] = lhs[cell.id()]*rhs;

    for(const Face &face: lhs.grid.faces())
        result.faces()[face.id()] = lhs.faces()[face.id()]*rhs;

    return result;
}

VectorFiniteVolumeField operator/(VectorFiniteVolumeField lhs, const ScalarFiniteVolumeField& rhs)
{
    lhs /= rhs;
    return lhs;
}
