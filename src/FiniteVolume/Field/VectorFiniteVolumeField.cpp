#include "VectorFiniteVolumeField.h"
#include "Exception.h"

template<>
VectorFiniteVolumeField& VectorFiniteVolumeField::operator=(const SparseVector& rhs)
{
    auto &self = *this;
    const size_t nActiveCells = grid.nActiveCells();

    for(int i = 0, end = nActiveCells; i < end; ++i)
        self[i].x = rhs[i];

    for(int i = nActiveCells, end = 2*nActiveCells; i < end; ++i)
        self[i - nActiveCells].y = rhs[i];

    return self;
}

//- Protected methods

template<>
void VectorFiniteVolumeField::setBoundaryRefValues(const Input &input)
{
    using namespace std;

    auto &self = *this;

    //- Check for a default patch

    auto defBoundary = input.boundaryInput().get_child_optional("Boundaries." + name + ".*");

    if(defBoundary)
    {
        for(int i = 0; i < grid.patches().size(); ++i)
            boundaryRefValues_.push_back(Vector2D(input.boundaryInput().get<std::string>("Boundaries." + name + ".*.value")));
    }
    else
    {
        //- Boundary condition association to patches is done by patch id
        for(const Patch& patch: grid.patches())
        {
            string root = "Boundaries." + name + "." + patch.name;
            boundaryRefValues_.push_back(Vector2D(input.boundaryInput().get<std::string>(root + ".value")));
        }
    }

    for(const Patch& patch: grid.patches())
    {
        for(const Face& face: patch.faces())
            self.faces()[face.id()] = boundaryRefValue(face.id());
    }
}

//- External functions

VectorFiniteVolumeField grad(const ScalarFiniteVolumeField &scalarField)
{
    VectorFiniteVolumeField gradField(scalarField.grid, "grad_" + scalarField.name);

    for(const Cell& cell: scalarField.grid.cells())
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
