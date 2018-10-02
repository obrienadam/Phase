#include "System/Exception.h"

#include "VectorFiniteVolumeField.h"

template<>
void VectorFiniteVolumeField::faceToCell(const CellGroup &cells)
{
    auto &self = *this;

    for (const Cell &cell: cells)
    {
        Vector2D sumSf(0., 0.), tmp(0., 0.);

        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.outwardNorm().abs();
            tmp += pointwise(self(nb.face()), sf);
            sumSf += sf;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.outwardNorm().abs();
            tmp += pointwise(self(bd.face()), sf);
            sumSf += sf;
        }

        self(cell) = Vector2D(tmp.x / sumSf.x, tmp.y / sumSf.y);
    }
}

template<>
void VectorFiniteVolumeField::faceToCell(const FiniteVolumeField<Scalar> &cellWeight,
                                         const FiniteVolumeField<Scalar> &faceWeight,
                                         const CellGroup &cells)
{
    auto &self = *this;

    for (const Cell &cell: cells)
    {
        Vector2D sumSf(0., 0.), tmp(0., 0.);

        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.outwardNorm().abs();
            tmp += pointwise(self(nb.face()), sf) / faceWeight(nb.face());
            sumSf += sf;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.outwardNorm().abs();
            tmp += pointwise(self(bd.face()), sf) / faceWeight(bd.face());
            sumSf += sf;
        }

        self(cell) = cellWeight(cell) * Vector2D(tmp.x / sumSf.x, tmp.y / sumSf.y);
    }
}

template<>
void VectorFiniteVolumeField::faceToCellAxisymmetric(const CellGroup &cells)
{
    auto &self = *this;

    for (const Cell &cell: cells)
    {
        Vector2D sumSf(0., 0.), tmp(0., 0.);

        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.polarOutwardNorm().abs();
            tmp += pointwise(self(nb.face()), sf);
            sumSf += sf;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.polarOutwardNorm().abs();
            tmp += pointwise(self(bd.face()), sf);
            sumSf += sf;
        }

        self(cell) = Vector2D(tmp.x / sumSf.x, tmp.y / sumSf.y);
    }
}

template<>
void VectorFiniteVolumeField::faceToCellAxisymmetric(const FiniteVolumeField<Scalar> &cw, const FiniteVolumeField<Scalar> &fw, const CellGroup &cells)
{
    auto &self = *this;

    for (const Cell &cell: cells)
    {
        Vector2D sumSf(0., 0.), tmp(0., 0.);

        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.polarOutwardNorm().abs();
            tmp += pointwise(self(nb.face()), sf) / fw(nb.face());
            sumSf += sf;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.polarOutwardNorm().abs();
            tmp += pointwise(self(bd.face()), sf) / fw(bd.face());
            sumSf += sf;
        }

        self(cell) = cw(cell) * Vector2D(tmp.x / sumSf.x, tmp.y / sumSf.y);
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

        for (const FaceGroup &patch: grid()->patches())
        {
            BoundaryType type = patchBoundaries_[patch.name()].first;
            patchBoundaries_[patch.name()] = std::make_pair(type, refVal);
        }
    }

    for (const FaceGroup &patch: grid()->patches())
    {
        valStr = input.boundaryInput().get<string>("Boundaries." + name_ + "." + patch.name() + ".value", "");

        if (valStr.empty())
            continue;

        Vector2D refVal(valStr);

        BoundaryType type = patchBoundaries_[patch.name()].first;
        patchBoundaries_[patch.name()] = std::make_pair(type, refVal);
    }

    auto &self = *this;
    for (const FaceGroup &patch: grid()->patches())
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

    for (const FaceGroup &patch: grid_->patches())
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
    VectorFiniteVolumeField result(lhs.grid(), lhs.name());

    for (const Cell &cell: lhs.grid()->cells())
        result[cell.id()] = lhs[cell.id()] * rhs;

    for (const Face &face: lhs.grid()->faces())
        result.faces()[face.id()] = lhs.faces()[face.id()] * rhs;

    return result;
}

VectorFiniteVolumeField operator/(VectorFiniteVolumeField lhs, const ScalarFiniteVolumeField &rhs)
{
    lhs /= rhs;
    return lhs;
}
