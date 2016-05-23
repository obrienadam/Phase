#include "Equation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

template<>
Equation<ScalarFiniteVolumeField>::Equation(ScalarFiniteVolumeField& field, const std::string& name)
    :
      name(name),
      spMat_(field.grid.nActiveCells(), field.grid.nActiveCells(), 5),
      field_(field)
{
    boundaries_ = SparseVector::Zero(field.grid.nActiveCells());
    sources_ = SparseVector::Zero(field.grid.nActiveCells());
}

template<>
Equation<VectorFiniteVolumeField>::Equation(VectorFiniteVolumeField& field, const std::string& name)
    :
      name(name),
      spMat_(2*field.grid.nActiveCells(), 2*field.grid.nActiveCells(), 5),
      field_(field)
{
    boundaries_ = SparseVector::Zero(2*field.grid.nActiveCells());
    sources_ = SparseVector::Zero(2*field.grid.nActiveCells());
}

template<>
Equation<ScalarFiniteVolumeField>& Equation<ScalarFiniteVolumeField>::operator +=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.fluidCells())
        sources_(cell.globalIndex()) += rhs[cell.id()];

    return *this;
}

template<>
Equation<ScalarFiniteVolumeField>& Equation<ScalarFiniteVolumeField>::operator -=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.fluidCells())
        sources_(cell.globalIndex()) -= rhs[cell.id()];

    return *this;
}

template<>
Equation<VectorFiniteVolumeField>& Equation<VectorFiniteVolumeField>::operator +=(const VectorFiniteVolumeField& rhs)
{
    const size_t nActiveCells = rhs.grid.nActiveCells();

    for(const Cell& cell: rhs.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        sources_(rowX) += rhs[cell.id()].x;
        sources_(rowY) += rhs[cell.id()].y;
    }

    return *this;
}

template<>
Equation<VectorFiniteVolumeField>& Equation<VectorFiniteVolumeField>::operator -=(const VectorFiniteVolumeField& rhs)
{
    const size_t nActiveCells = rhs.grid.nActiveCells();

    for(const Cell& cell: rhs.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        sources_(rowX) -= rhs[cell.id()].x;
        sources_(rowY) -= rhs[cell.id()].y;
    }

    return *this;
}

template<>
void Equation<ScalarFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    for(const Cell& cell: field_.grid.fluidCells()) // Should the whole matrix be relaxed??
    {
        size_t idx = cell.globalIndex();

        spMat_.coeffRef(idx, idx) /= relaxationFactor;
        boundaries_(idx) += (1. - relaxationFactor)*spMat_.coeff(idx, idx)*field_[cell.id()];
    }
}

template<>
void Equation<VectorFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    const size_t nActiveCells = field_.grid.nActiveCells();

    for(const Cell& cell: field_.grid.fluidCells()) // Should the whole matrix be relaxed??
    {
        size_t idxX = cell.globalIndex();
        size_t idxY = idxX + nActiveCells;

        spMat_.coeffRef(idxX, idxX) /= relaxationFactor;
        boundaries_(idxX) += (1. - relaxationFactor)*spMat_.coeff(idxX, idxX)*field_[cell.id()].x;

        spMat_.coeffRef(idxY, idxY) /= relaxationFactor;
        boundaries_(idxY) += (1. - relaxationFactor)*spMat_.coeff(idxY, idxY)*field_[cell.id()].y;
    }
}

//- External functions

namespace fv
{
Equation<ScalarFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& gamma, ScalarFiniteVolumeField& field)
{
    const size_t nCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(5*nCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t row = cell.globalIndex();
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            size_t col = nb.cell().globalIndex();
            Scalar coeff = gamma.faces()[nb.face().id()]*dot(nb.rCellVec(), nb.outwardNorm())/dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar coeff = gamma.faces()[bd.face().id()]*dot(bd.rFaceVec(), bd.outwardNorm())/dot(bd.rFaceVec(), bd.rFaceVec());

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                eqn.boundaries()(row) -= coeff*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                eqn.boundaries()(row) -= bd.outwardNorm().mag()/coeff*field.boundaryRefValue(bd.face().id());
                break;

            default:
                throw Exception("fv", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& gamma, VectorFiniteVolumeField& field)
{
    const size_t nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<VectorFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);

    entries.reserve(10*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            size_t colX = nb.cell().globalIndex();
            size_t colY = colX + nActiveCells;

            Scalar coeff = gamma.faces()[nb.face().id()]*dot(nb.rCellVec(), nb.outwardNorm())/dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;

            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, colX, coeff));
            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, colY, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar coeff = gamma.faces()[bd.face().id()]*dot(bd.rFaceVec(), bd.outwardNorm())/dot(bd.rFaceVec(), bd.rFaceVec());
            Vector2D source;

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                eqn.boundaries()(rowX) -= coeff*field.faces()[bd.face().id()].x; // THIS IS THE MOTHER FUCKING PROBLEM O_O
                eqn.boundaries()(rowY) -= coeff*field.faces()[bd.face().id()].y;
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                source = bd.outwardNorm().mag()/coeff*field.boundaryRefValue(bd.face().id());
                eqn.boundaries()(rowX) -= source.x;
                eqn.boundaries()(rowY) -= source.y;
                break;

            default:
                throw Exception("fv", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, centralCoeff));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField& u, ScalarFiniteVolumeField& field)
{
    const size_t nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(5*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t row = cell.globalIndex();
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            size_t col = nb.cell().globalIndex();

            Scalar faceFlux = dot(u.faces()[nb.face().id()], nb.outwardNorm());

            Scalar coeff = std::min(faceFlux, 0.);
            centralCoeff += std::max(faceFlux, 0.);

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar faceFlux = dot(u.faces()[bd.face().id()], bd.outwardNorm());

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= faceFlux*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += faceFlux;
                eqn.boundaries()(row) -= bd.outwardNorm().mag()*field.boundaryRefValue(bd.face().id());
                break;

            default:
                throw Exception("fv", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(row, row, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> div(const VectorFiniteVolumeField& u, VectorFiniteVolumeField& field)
{
    const size_t nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);

    entries.reserve(5*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            size_t colX = nb.cell().globalIndex();
            size_t colY = colX + nActiveCells;

            Scalar faceFlux = dot(u.faces()[nb.face().id()], nb.outwardNorm());

            Scalar coeff = std::min(faceFlux, 0.);
            centralCoeff += std::max(faceFlux, 0.);

            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, colX, coeff));
            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, colY, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar faceFlux = dot(u.faces()[bd.face().id()], bd.outwardNorm());

            switch(field.boundaryType(bd.face().id()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.boundaries()(rowX) -= faceFlux*field.faces()[bd.face().id()].x;
                eqn.boundaries()(rowY) -= faceFlux*field.faces()[bd.face().id()].y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += faceFlux;
                eqn.boundaries()(rowX) -= bd.outwardNorm().mag()*field.boundaryRefValue(bd.face().id()).x;
                eqn.boundaries()(rowY) -= bd.outwardNorm().mag()*field.boundaryRefValue(bd.face().id()).y;
                break;

            default:
                throw Exception("fv", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, centralCoeff));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<ScalarFiniteVolumeField> ddt(const ScalarFiniteVolumeField& a, ScalarFiniteVolumeField& field, Scalar timeStep, Scalar prevTimeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();
    const Field<Scalar> &prevField = field.prev(0);

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t row = cell.globalIndex();

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, a[cell.id()]*cell.volume()/timeStep));
        eqn.boundaries()(row) += a[cell.id()]*cell.volume()*prevField[cell.id()]/timeStep;
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<ScalarFiniteVolumeField> ddt(ScalarFiniteVolumeField& field, Scalar timeStep, Scalar prevTimeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();
    const Field<Scalar> &prevField = field.prev(0);

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t row = cell.globalIndex();

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, cell.volume()/timeStep));
        eqn.boundaries()(row) += cell.volume()*prevField[cell.id()]/timeStep;
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> ddt(VectorFiniteVolumeField& field, Scalar timeStep, Scalar prevTimeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();
    const FiniteVolumeField<Vector2D> &prevField = field.prev(0);

    std::vector<Equation<VectorFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);

    entries.reserve(2*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, cell.volume()/timeStep));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, cell.volume()/timeStep));

        eqn.boundaries()(rowX) += cell.volume()*prevField[cell.id()].x/timeStep;
        eqn.boundaries()(rowY) += cell.volume()*prevField[cell.id()].y/timeStep;
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> ddt(const ScalarFiniteVolumeField& a, VectorFiniteVolumeField& field, Scalar timeStep, Scalar prevTimeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();
    const FiniteVolumeField<Vector2D> &prevField = field.prev(0);

    std::vector<Equation<VectorFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);

    entries.reserve(2*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, a[cell.id()]*cell.volume()/timeStep));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, a[cell.id()]*cell.volume()/timeStep));

        eqn.boundaries()(rowX) += a[cell.id()]*cell.volume()*prevField[cell.id()].x/timeStep;
        eqn.boundaries()(rowY) += a[cell.id()]*cell.volume()*prevField[cell.id()].y/timeStep;
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

VectorFiniteVolumeField grad(const ScalarFiniteVolumeField &field)
{
    VectorFiniteVolumeField gradField(field.grid, "grad_" + field.name);

    for(const Cell& cell: gradField.grid.fluidCells())
    {
        Vector2D &gradVec = gradField[cell.id()];

        for(const InteriorLink& nb: cell.neighbours())
            gradVec += field.faces()[nb.face().id()]*nb.outwardNorm();

        for(const BoundaryLink& bd: cell.boundaries())
            gradVec += field.faces()[bd.face().id()]*bd.outwardNorm();
    }

    return gradField;
}

VectorFiniteVolumeField source(VectorFiniteVolumeField field)
{
    for(const Cell &cell: field.grid.fluidCells())
        field[cell.id()] *= cell.volume();

    return field;
}

}
