#include "VectorFiniteVolumeField.h"
#include "Exception.h"

template<>
VectorFiniteVolumeField& VectorFiniteVolumeField::operator=(const SparseVector& rhs)
{
    auto &self = *this;
    const size_t nActiveCells = grid.nActiveCells();

    for(const Cell &cell: grid.activeCells())
    {
        self[cell.id()].x = rhs[cell.globalIndex()];
        self[cell.id()].y = rhs[cell.globalIndex() + nActiveCells];
    }

    return self;
}

//- Protected methods

template<>
void VectorFiniteVolumeField::setBoundaryRefValues(const Input &input)
{
    using namespace std;

    string valStr = input.boundaryInput().get<string>("Boundaries." + name + ".*.value", "");

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
        valStr = input.boundaryInput().get<string>("Boundaries." + name + "." + entry.first + ".value", "");

        if(valStr.empty())
            continue;

        Vector2D refVal(valStr);

        BoundaryType type = patchBoundaries_[entry.second.id()].first;
        patchBoundaries_[entry.second.id()] = std::make_pair(type, refVal);
    }

    auto &self = *this;

    for(const Face &face: grid.boundaryFaces())
        self.faces()[face.id()] = boundaryRefValue(face.id());
}

//- External functions

VectorFiniteVolumeField grad(const ScalarFiniteVolumeField &scalarField)
{
    VectorFiniteVolumeField gradField(scalarField.grid, "grad_" + scalarField.name);

    for(const Cell& cell: scalarField.grid.activeCells())
    {
        Vector2D &gradPhi = gradField[cell.id()];

        for(const InteriorLink& nb: cell.neighbours())
            gradPhi += scalarField.faces()[nb.face().id()]*nb.outwardNorm();

        for(const BoundaryLink& bd: cell.boundaries())
            gradPhi += scalarField.faces()[bd.face().id()]*bd.outwardNorm();

        gradPhi /= cell.volume();
    }

    return gradField;
}

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
    VectorFiniteVolumeField result(lhs.grid, lhs.name);

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
