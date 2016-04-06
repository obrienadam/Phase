#include "Equation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

template<>
Equation<ScalarFiniteVolumeField>::Equation(ScalarFiniteVolumeField& field)
    :
      spMat_(field.grid.nActiveCells(), field.grid.nActiveCells(), 5),
      b_(field.grid.nActiveCells()),
      field_(field)
{
    b_.setZero();
}

template<>
Equation<VectorFiniteVolumeField>::Equation(VectorFiniteVolumeField& field)
    :
      spMat_(2*field.grid.nActiveCells(), 2*field.grid.nActiveCells(), 5),
      b_(2*field.grid.nActiveCells()),
      field_(field)
{
    b_.setZero();
}

template<>
Equation<ScalarFiniteVolumeField>& Equation<ScalarFiniteVolumeField>::operator +=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.cells)
    {
        if(!cell.isActive())
            continue;

        b_[cell.globalIndex()] += rhs[cell.id()];
    }

    return *this;
}

template<>
Equation<ScalarFiniteVolumeField>& Equation<ScalarFiniteVolumeField>::operator -=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.cells)
    {
        if(!cell.isActive())
            continue;

        b_[cell.globalIndex()] -= rhs[cell.id()];
    }

    return *this;
}

template<>
Equation<VectorFiniteVolumeField>& Equation<VectorFiniteVolumeField>::operator +=(const VectorFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.cells)
    {
        if(!cell.isActive())
            continue;

        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + rhs.grid.nActiveCells();

        b_[rowX] += rhs[cell.id()].x;
        b_[rowY] += rhs[cell.id()].y;
    }

    return *this;
}

template<>
Equation<VectorFiniteVolumeField>& Equation<VectorFiniteVolumeField>::operator -=(const VectorFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.cells)
    {
        if(!cell.isActive())
            continue;

        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + rhs.grid.nActiveCells();

        b_[rowX] -= rhs[cell.id()].x;
        b_[rowY] -= rhs[cell.id()].y;
    }

    return *this;
}

template<>
void Equation<ScalarFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    for(const Cell& cell: field_.grid.cells)
    {
        if(!cell.isActive())
            continue;

        size_t idx = cell.globalIndex();

        spMat_.coeffRef(idx, idx) /= relaxationFactor;
        b_(idx) += (1. - relaxationFactor)*spMat_.coeff(idx, idx)*field_[cell.id()];
    }
}

template<>
void Equation<VectorFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    const size_t nActiveCells = field_.grid.nActiveCells();

    for(const Cell& cell: field_.grid.cells)
    {
        if(!cell.isActive())
            continue;

        size_t idxX = cell.globalIndex();
        size_t idxY = idxX + nActiveCells;

        spMat_.coeffRef(idxX, idxX) /= relaxationFactor;
        b_(idxX) += (1. - relaxationFactor)*spMat_.coeff(idxX, idxX)*field_[cell.id()].x;

        spMat_.coeffRef(idxY, idxY) /= relaxationFactor;
        b_(idxY) += (1. - relaxationFactor)*spMat_.coeff(idxY, idxY)*field_[cell.id()].y;
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

    for(const Cell& cell: field.grid.cells)
    {
        if(!cell.isActive())
            continue;

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
                eqn.rhs()[row] = -coeff*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                break;

            default:
                throw Exception("", "laplacian", "unrecognized or unspecified boundary type.");
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

    for(const Cell& cell: field.grid.cells)
    {
        if(!cell.isActive())
            continue;

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

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                eqn.rhs()(rowX) = -coeff*field.faces()[bd.face().id()].x;
                eqn.rhs()(rowY) = -coeff*field.faces()[bd.face().id()].y;
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                break;

            default:
                throw Exception("", "laplacian", "unrecognized or unspecified boundary type.");
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

    for(const Cell& cell: field.grid.cells)
    {
        if(!cell.isActive())
            continue;

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
            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.rhs()[row] = -dot(u.faces()[bd.face().id()], bd.outwardNorm())*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                break;

            default:
                throw Exception("", "div", "unrecognized or unspecified boundary type.");
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

    for(const Cell& cell: field.grid.cells)
    {
        if(!cell.isActive())
            continue;

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
            switch(field.boundaryType(bd.face().id()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.rhs()[rowX] = -dot(u.faces()[bd.face().id()], bd.outwardNorm())*field.faces()[bd.face().id()].x;
                eqn.rhs()[rowY] = -dot(u.faces()[bd.face().id()], bd.outwardNorm())*field.faces()[bd.face().id()].y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                break;

            default:
                throw Exception("", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, centralCoeff));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<ScalarFiniteVolumeField> ddt(const ScalarFiniteVolumeField& a, ScalarFiniteVolumeField& field, Scalar timeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(nActiveCells);

    for(const Cell& cell: field.grid.cells)
    {
        if(!cell.isActive())
            continue;

        size_t row = cell.globalIndex();

        Scalar coeff = a[cell.id()]*cell.volume()/timeStep;

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, coeff));
        eqn.rhs()[row] = coeff*field[cell.id()];
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> ddt(const ScalarFiniteVolumeField& a, VectorFiniteVolumeField& field, Scalar timeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<VectorFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);

    entries.reserve(2*nActiveCells);

    for(const Cell& cell: field.grid.cells)
    {
        if(!cell.isActive())
            continue;

        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        Scalar coeff = a[cell.id()]*cell.volume()/timeStep;

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, coeff));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, coeff));

        eqn.rhs()[rowX] = coeff*field[cell.id()].x;
        eqn.rhs()[rowY] = coeff*field[cell.id()].y;
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

VectorFiniteVolumeField grad(const ScalarFiniteVolumeField &field)
{
    VectorFiniteVolumeField gradField(field.grid, "grad_" + field.name);

    for(const Cell& cell: gradField.grid.cells)
    {
        if(!cell.isActive())
            continue;

        Vector2D &gradVec = gradField[cell.id()];

        for(const InteriorLink& nb: cell.neighbours())
            gradVec += field.faces()[nb.face().id()]*nb.outwardNorm();

        for(const BoundaryLink& bd: cell.boundaries())
            gradVec += field.faces()[bd.face().id()]*bd.outwardNorm();
    }

    return gradField;
}
}
