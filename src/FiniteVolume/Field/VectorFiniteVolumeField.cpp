#include "VectorFiniteVolumeField.h"
#include "FiniteVolumeGrid2D.h"
#include "Exception.h"

VectorFiniteVolumeField::VectorFiniteVolumeField(const FiniteVolumeGrid2D &grid, const std::string &name)
    :
      Field<Vector2D>::Field(grid.cells.size(), Vector2D(), name),
      faces_(grid.faces.size(), Vector2D()),
      grid(grid)
{

}

VectorFiniteVolumeField::VectorFiniteVolumeField(const Input &input, const FiniteVolumeGrid2D &grid, const std::string &name)
    :
      VectorFiniteVolumeField(grid, name)
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
            throw Exception("VectorFiniteVolumeField", "VectorFiniteVolumeField", "unrecognized boundary type \"" + type + "\".");
        }

        std::string refValStr = input.boundaryInput().get<std::string>(root + ".value");
        boundaryRefValues_.push_back(Vector2D(refValStr));

        for(const Face& face: patch.faces())
        {
            self.faces()[face.id()] = boundaryRefValues_.back();
        }
    }
}

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

VectorFiniteVolumeField& VectorFiniteVolumeField::operator =(const VectorFiniteVolumeField& rhs)
{
    if(this == &rhs)
        return *this;

    assert(&grid == &rhs.grid);

    boundaryTypes_ = rhs.boundaryTypes_;
    boundaryRefValues_ = rhs.boundaryRefValues_;
    faces_ = rhs.faces_;
    Field<Vector2D>::operator =(rhs);

    return *this;
}

VectorFiniteVolumeField& VectorFiniteVolumeField::operator *=(const ScalarFiniteVolumeField& rhs)
{
    auto &self = *this;

    for(int i = 0, end = self.size(); i < end; ++i)
        self[i] *= rhs[i];

    for(int i = 0, end = self.faces().size(); i < end; ++i)
        self.faces()[i] *= rhs.faces()[i];

    return self;
}

VectorFiniteVolumeField::BoundaryType VectorFiniteVolumeField::boundaryType(size_t faceId) const
{
    if(boundaryTypes_.size() == 0)
        return VectorFiniteVolumeField::NORMAL_GRADIENT;

    return boundaryTypes_[grid.faces[faceId].patch().id()];
}

//- External functions

VectorFiniteVolumeField grad(const ScalarFiniteVolumeField &scalarField)
{
    VectorFiniteVolumeField gradField(scalarField.grid, "grad_" + scalarField.name);

    for(const Cell& cell: scalarField.grid.cells)
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

void interpolateFaces(VectorFiniteVolumeField& field)
{
    for(const Face& face: field.grid.faces)
    {
        if(face.isBoundary())
            continue;

        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        Scalar alpha = rCell.volume()/(lCell.volume() + rCell.volume());
        field.faces()[face.id()] = field[lCell.id()]*alpha + field[rCell.id()]*(1. - alpha);
    }
}

VectorFiniteVolumeField operator+(VectorFiniteVolumeField lhs, const VectorFiniteVolumeField& rhs)
{
    lhs += rhs;
    return lhs;
}

VectorFiniteVolumeField operator-(VectorFiniteVolumeField lhs, const VectorFiniteVolumeField& rhs)
{
    lhs -= rhs;
    return lhs;
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
