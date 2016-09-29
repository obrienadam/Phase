#include "TotalVariationDiminishing.h"
#include "GradientEvaluation.h"

namespace tvd
{

Scalar superbee(Scalar r)
{
    return std::max(0., std::max(std::min(1., 2.*r), std::min(2., r)));
}

Scalar minmod(Scalar r)
{
    return std::max(0., std::min(1., r));
}

Scalar osher(Scalar r)
{
    return std::max(0., std::min(2., r));
}

Scalar muscl(Scalar r)
{
    return (r + fabs(r))/(1. + fabs(r));
}

Scalar charm(Scalar r)
{
    return r > 0. ? r*(3.*r + 1.)/(r*r + 2.*r + 1.) : 0.;
}

Equation<VectorFiniteVolumeField> div(const VectorFiniteVolumeField &u, VectorFiniteVolumeField &field, Limiter limiter)
{
    typedef Equation<VectorFiniteVolumeField>::Triplet Triplet;

    const Size nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);
    TensorFiniteVolumeField gradField(field.grid, "gradField");
    fv::computeGradient(fv::FACE_TO_CELL, field, gradField);

    entries.reserve(10*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        const Index rowX = cell.globalIndex();
        const Index rowY = rowX + nActiveCells;
        Scalar centralCoeffX = 0., centralCoeffY = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Index colX = nb.cell().globalIndex();
            const Index colY = colX + nActiveCells;
            const Vector2D& rc = nb.rCellVec();

            Vector2D r;

            Scalar flux = dot(u(nb.face()), nb.outwardNorm());

            if(flux > 0)
            {
                Vector2D tmp = 2*dot(gradField(nb.cell()), rc);
                r.x = tmp.x/(field(cell).x - field(nb.cell()).x) - 1;
                r.y = tmp.y/(field(cell).y - field(nb.cell()).y) - 1;
            }
            else
            {
                Vector2D tmp = 2*dot(gradField(cell), rc);
                r.x = tmp.x/(field(nb.cell()).x - field(cell).x) - 1;
                r.y = tmp.y/(field(nb.cell()).y - field(cell).y) - 1;
            }

            r.x = isnan(r.x) ? 0 : r.x;
            r.y = isnan(r.y) ? 0 : r.y;

            Vector2D phi;

            switch(limiter)
            {
            case SUPERBEE:
                phi = Vector2D(muscl(r.x), muscl(r.y));
                break;

            case MINMOD:
                phi = Vector2D(minmod(r.x), muscl(r.y));
                break;

            case OSHER:
                phi = Vector2D(osher(r.x), osher(r.y));
                break;

            case MUSCL:
                phi = Vector2D(muscl(r.x), muscl(r.y));
                break;

            case CHARM:
                phi = Vector2D(charm(r.x), charm(r.y));
                break;
            }

            phi.x = isnan(phi.x) ? 0 : phi.x;
            phi.y = isnan(phi.y) ? 0 : phi.y;

            const Scalar faceFlux = dot(u(nb.face()), nb.outwardNorm());

            const Scalar coeffX = faceFlux*phi.x/2.;
            const Scalar coeffY = faceFlux*phi.y/2.;
            centralCoeffX += faceFlux*(1. - phi.x/2.);
            centralCoeffY += faceFlux*(1. - phi.y/2.);

            entries.push_back(Triplet(rowX, colX, coeffX));
            entries.push_back(Triplet(rowY, colY, coeffY));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());

            switch(field.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.boundaries()(rowX) -= faceFlux*field(bd.face()).x;
                eqn.boundaries()(rowY) -= faceFlux*field(bd.face()).y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeffX += faceFlux;
                centralCoeffY += faceFlux;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                Scalar boundaryCoeff = dot(bd.rFaceVec(), bd.outwardNorm().unitVec())/dot(bd.rFaceVec(), bd.rFaceVec());
                boundaryCoeff = 1./(1./boundaryCoeff + 1.);
                centralCoeffX += faceFlux*boundaryCoeff;
                centralCoeffY += faceFlux*boundaryCoeff;
            }
                break;

            default:
                throw Exception("tvd", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, centralCoeffX));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, centralCoeffY));
    }

    eqn.assemble(entries);
    return eqn;
}

}
