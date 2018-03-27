#include "ScalarGradient.h"

ScalarGradient::ScalarGradient(const ScalarFiniteVolumeField &phi)
        :
        VectorFiniteVolumeField(phi.grid(),
                                "grad" + phi.name(),
                                Vector2D(0., 0.),
                                true, false),
        phi_(phi)
{

}

void ScalarGradient::computeFaces()
{
    VectorFiniteVolumeField &gradPhi = *this;

    for (const Face &face: grid_->interiorFaces())
    {
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        gradPhi(face) = (phi_(face.rCell()) - phi_(face.lCell())) * rc / dot(rc, rc);
    }

    for (const Face &face: grid_->boundaryFaces())
    {
        Vector2D rf = face.centroid() - face.lCell().centroid();
        gradPhi(face) = (phi_(face) - phi_(face.lCell())) * rf / dot(rf, rf);
    }
}

void ScalarGradient::compute(const CellGroup &group, Method method)
{
    computeFaces();
    VectorFiniteVolumeField &gradPhi = *this;

    //std::fill(gradPhi.begin(), gradPhi.end(), Vector2D(0., 0.));

    switch (method)
    {
        case FACE_TO_CELL:
            for (const Cell &cell: group)
            {
                Vector2D sum(0., 0.), tmp(0., 0.);

                for (const InteriorLink &nb: cell.neighbours())
                {
                    Vector2D sf = nb.outwardNorm().abs();
                    tmp += pointwise(gradPhi(nb.face()), sf);
                    sum += sf;
                }

                for (const InteriorLink &bd: cell.neighbours())
                {
                    Vector2D sf = bd.outwardNorm().abs();
                    tmp += pointwise(gradPhi(bd.face()), sf);
                    sum += sf;
                }

                gradPhi(cell) = Vector2D(tmp.x / sum.x, tmp.y / sum.y);
            }
            break;
        case GREEN_GAUSS_CELL:
            for (const Cell &cell: group)
            {
                for (const InteriorLink &nb: cell.neighbours())
                {
                    Scalar g = nb.distanceWeight();
                    Scalar phiF = g * phi_(cell) + (1. - g) * phi_(nb.cell());
                    gradPhi(cell) += phiF * nb.outwardNorm();
                }

                for (const BoundaryLink &bd: cell.boundaries())
                    gradPhi(cell) += phi_(bd.face()) * bd.outwardNorm();

                gradPhi(cell) /= cell.volume();
            }
            break;
        case GREEN_GAUSS_NODE:
            for (const Cell &cell: group)
            {
                for (const InteriorLink &nb: cell.neighbours())
                {
                    auto lNodeWeights = nb.face().lNode().distanceWeights();
                    auto rNodeWeights = nb.face().rNode().distanceWeights();
                    Scalar phiLN = 0, phiRN = 0;
                    int i = 0;
                    for (const Cell &cell: nb.face().lNode().cells())
                        phiLN += lNodeWeights[i++] * phi_(cell);

                    i = 0;
                    for (const Cell &cell: nb.face().rNode().cells())
                        phiRN += rNodeWeights[i++] * phi_(cell);

                    Scalar phiF = (phiLN + phiRN) / 2.;
                    gradPhi(cell) += phiF * nb.outwardNorm();
                }

                for (const BoundaryLink &bd: cell.boundaries())
                {
                    auto lNodeWeights = bd.face().lNode().distanceWeights();
                    auto rNodeWeights = bd.face().rNode().distanceWeights();
                    Scalar phiLN = 0, phiRN = 0;
                    int i = 0;
                    for (const Cell &cell: bd.face().lNode().cells())
                        phiLN += lNodeWeights[i++] * phi_(cell);

                    i = 0;
                    for (const Cell &cell: bd.face().rNode().cells())
                        phiRN += rNodeWeights[i++] * phi_(cell);

                    Scalar phiF = (phiLN + phiRN) / 2.;
                    gradPhi(cell) += phiF * bd.outwardNorm();
                }

                gradPhi(cell) /= cell.volume();
            }
            break;
    }
}

void ScalarGradient::compute(Method method)
{
    compute(phi_.cells(), method);
}

void ScalarGradient::computeAxisymmetric(const CellGroup &cells,
                                         Method method)
{
    computeFaces();
    std::fill(this->begin(), this->end(), Vector2D(0., 0.));
    VectorFiniteVolumeField &gradPhi = *this;

    switch (method)
    {
        case FACE_TO_CELL:
            for (const Cell &cell: cells)
            {
                Vector2D sum(0., 0.), tmp(0., 0.);

                for (const InteriorLink &nb: cell.neighbours())
                {
                    Vector2D sf = nb.face().polarOutwardNorm(cell.centroid()).abs();
                    tmp += pointwise(gradPhi(nb.face()), sf);
                    sum += sf;
                }

                for (const InteriorLink &bd: cell.neighbours())
                {
                    Vector2D sf = bd.face().polarOutwardNorm(cell.centroid()).abs();
                    tmp += pointwise(gradPhi(bd.face()), sf);
                    sum += sf;
                }

                gradPhi(cell) = Vector2D(tmp.x / sum.x, tmp.y / sum.y);
            }
            break;
        case GREEN_GAUSS_CELL:
//            for (const Cell &cell: cells)
//            {
//                for (const InteriorLink &nb: cell.neighbours())
//                {
//                    Scalar g = nb.distanceWeight();
//                    Scalar phiF = g * phi_(cell) + (1. - g) * phi_(nb.cell());
//                    gradPhi(cell) += phiF * nb.face().polarOutwardNorm(cell.centroid(), axis, axisPoint);
//                }
//
//                for (const BoundaryLink &bd: cell.boundaries())
//                    gradPhi(cell) += phi_(bd.face()) * bd.face().polarOutwardNorm(cell.centroid(), axis, axisPoint);
//
//                gradPhi(cell) /= cell.polarVolume();
//            }
            break;
    }
}