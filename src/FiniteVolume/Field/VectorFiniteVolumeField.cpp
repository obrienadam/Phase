#include "VectorFiniteVolumeField.h"
#include "Exception.h"

template<>
VectorFiniteVolumeField &VectorFiniteVolumeField::operator=(const Vector &rhs)
{
    auto &self = *this;
    const size_t nActiveCells = grid().nLocalActiveCells();

    for (const Cell &cell: grid().localActiveCells())
    {
        self[cell.id()].x = rhs[cell.index(0)];
        self[cell.id()].y = rhs[cell.index(0) + nActiveCells];
    }

    return self;
}

template<>
Vector VectorFiniteVolumeField::vectorize() const
{
    const auto &self = *this;
    const Size nActiveCells = grid().nLocalActiveCells();
    Vector vec(2 * nActiveCells, 0.);

    for (const Cell &cell: grid().localActiveCells())
    {
        vec[cell.index(0)] = self(cell).x;
        vec[nActiveCells + cell.index(0)] = self(cell).y;
    }

    return vec;
}

template<>
void VectorFiniteVolumeField::faceToCell(const FiniteVolumeField<Scalar> &cellWeight,
                                         const FiniteVolumeField<Scalar> &faceWeight)
{
    auto &self = *this;

    for(const Cell& cell: grid_->cells())
    {
        Vector2D sumSf(0., 0.), tmp(0., 0.);

        for(const InteriorLink& nb: cell.neighbours())
        {
            Vector2D sf = nb.outwardNorm().abs();
            tmp += pointwise(self(nb.face()), sf) / faceWeight(nb.face());
            sumSf += sf;
        }

        for(const BoundaryLink& bd: cell.boundaries())
        {
            Vector2D sf = bd.outwardNorm().abs();
            tmp += pointwise(self(bd.face()), sf) / faceWeight(bd.face());
            sumSf += sf;
        }

        self(cell) = cellWeight(cell) * Vector2D(tmp.x / sumSf.x, tmp.y / sumSf.y);
    }
}

//- Protected methods

template<>
void VectorFiniteVolumeField::setBoundaryRefValues(const Input &input)
{
    using namespace std;

    string valStr = input.boundaryInput().get<string>("Boundaries." + name_ + ".*.value", "");

    if (!valStr.empty())
    {
        Vector2D refVal = Vector2D(valStr);

        for (const Patch &patch: grid().patches())
        {
            BoundaryType type = patchBoundaries_[patch.id()].first;
            patchBoundaries_[patch.id()] = std::make_pair(type, refVal);
        }
    }

    for (const Patch &patch: grid().patches())
    {
        valStr = input.boundaryInput().get<string>("Boundaries." + name_ + "." + patch.name() + ".value", "");

        if (valStr.empty())
            continue;

        Vector2D refVal(valStr);

        BoundaryType type = patchBoundaries_[patch.id()].first;
        patchBoundaries_[patch.id()] = std::make_pair(type, refVal);
    }

    auto &self = *this;
    for (const Patch &patch: grid().patches())
    {
        Vector2D bRefVal = boundaryRefValue(patch);

        for (const Face &face: patch)
            self(face) = bRefVal;
    }
}

template<>
void VectorFiniteVolumeField::setBoundaryFaces()
{
    auto &self = *this;

    for (const Patch &patch: grid_->patches())
        switch (boundaryType(patch))
        {
            case FIXED:
                break;
            case NORMAL_GRADIENT:
                for (const Face &face: patch)
                    self(face) = self(face.lCell());
                break;
            case SYMMETRY:
                for (const Face &face: patch)
                    self(face) = self(face.lCell())
                                 - dot(self(face.lCell()), face.norm()) * face.norm() / face.norm().magSqr();
                break;
            default:
                throw Exception("VectorFiniteVolumeField", "setBoundaryFaces", "unrecognized boundary type.");
        }
}

//- External functions

VectorFiniteVolumeField operator*(const ScalarFiniteVolumeField &lhs, VectorFiniteVolumeField rhs)
{
    rhs *= lhs;
    return rhs;
}

VectorFiniteVolumeField operator*(VectorFiniteVolumeField lhs, const ScalarFiniteVolumeField &rhs)
{
    lhs *= rhs;
    return lhs;
}

VectorFiniteVolumeField operator*(const ScalarFiniteVolumeField &lhs, const Vector2D &rhs)
{
    VectorFiniteVolumeField result(lhs.gridPtr(), lhs.name());

    for (const Cell &cell: lhs.grid().cells())
        result[cell.id()] = lhs[cell.id()] * rhs;

    for (const Face &face: lhs.grid().faces())
        result.faces()[face.id()] = lhs.faces()[face.id()] * rhs;

    return result;
}

VectorFiniteVolumeField operator/(VectorFiniteVolumeField lhs, const ScalarFiniteVolumeField &rhs)
{
    lhs /= rhs;
    return lhs;
}
