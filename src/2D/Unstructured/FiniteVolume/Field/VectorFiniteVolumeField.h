#ifndef PHASE_VECTOR_FINITE_VOLUME_FIELD_H
#define PHASE_VECTOR_FINITE_VOLUME_FIELD_H

#include "Geometry/Vector2D.h"
#include "ScalarFiniteVolumeField.h"

typedef FiniteVolumeField<Vector2D> VectorFiniteVolumeField;

//- Specializations

template<>
template<class UnaryPredicate>
void VectorFiniteVolumeField::faceToCell(const FiniteVolumeField<Scalar> &cellWeight,
                                         const FiniteVolumeField<Scalar> &faceWeight,
                                         const CellGroup &cells,
                                         const UnaryPredicate &p)
{
    auto &self = *this;

    for (const Cell &cell: cells)
    {
        if(!p(cell))
            continue;

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
void VectorFiniteVolumeField::faceToCell(const FiniteVolumeField<Scalar> &cellWeight,
                                         const FiniteVolumeField<Scalar> &faceWeight,
                                         const CellGroup &cells);

template<>
void VectorFiniteVolumeField::setBoundaryRefValues(const Input &input);

template<>
void VectorFiniteVolumeField::setBoundaryFaces();

//- External
VectorFiniteVolumeField operator*(const ScalarFiniteVolumeField &lhs, VectorFiniteVolumeField rhs);

VectorFiniteVolumeField operator*(VectorFiniteVolumeField lhs, const ScalarFiniteVolumeField &rhs);

VectorFiniteVolumeField operator*(const ScalarFiniteVolumeField &lhs, const Vector2D &rhs);

VectorFiniteVolumeField operator/(VectorFiniteVolumeField lhs, const ScalarFiniteVolumeField &rhs);

#endif
