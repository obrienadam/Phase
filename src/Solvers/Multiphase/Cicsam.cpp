#include "Cicsam.h"

namespace cicsam
{

Scalar hc(Scalar gammaTilde, Scalar coD) // Maintains very sharp interface but funny things sometimes happen
{
    return gammaTilde >= 0 && gammaTilde <= 1 ? std::min(1., gammaTilde/coD) : gammaTilde;
}


Scalar uq(Scalar gammaTilde, Scalar coD) // Fairly sharp interface but fewer funny things
{
    return gammaTilde >= 0 && gammaTilde <= 1 ? std::min((8.*coD*gammaTilde + (1. - coD)*(6.*gammaTilde + 3.))/8., hc(gammaTilde, coD)) : gammaTilde;
}

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field, Scalar timeStep, Type type)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);
    VectorFiniteVolumeField gradField = grad(field);

    entries.reserve(5*field.grid.nActiveCells());

    for(const Cell &cell: field.grid.fluidCells())
    {
        Scalar centralCoeff = 0.;
        const size_t row = cell.globalIndex();

        for(const InteriorLink &nb: cell.neighbours())
        {
            //- Determine the donor and acceptor cell

            const Scalar flux = dot(u.faces()[nb.face().id()], nb.outwardNorm());
            const Cell &donor = flux > 0. ? cell : nb.cell();
            const Cell &acceptor = flux > 0. ? nb.cell() : cell;

            const Scalar gammaD = field[donor.id()];
            const Scalar gammaA = field[acceptor.id()];
            const Scalar gammaU = std::max(std::min(gammaA + 2*dot(donor.centroid() - acceptor.centroid(), gradField[donor.id()]), 1.), 0.);

            const Scalar coD = sqrt(u.faces()[nb.face().id()].magSqr()*timeStep*timeStep/nb.rCellVec().magSqr());

            const Scalar gammaTilde = (gammaD - gammaU)/(gammaA - gammaU);

            Scalar gammaTildeF;

            switch(type)
            {
            case HC:
                gammaTildeF = hc(gammaTilde, coD);
                break;
            case UQ:
                gammaTildeF = uq(gammaTilde, coD);
                break;
            }

            Scalar betaFace = (gammaTildeF - gammaTilde)/(1. - gammaTilde);

            if(isnan(betaFace)) // Upwind if a valid value is not found
                betaFace = &cell == &donor ? 0. : 1.;

            betaFace = std::max(std::min(betaFace, 1.), 0.);

            const size_t col = nb.cell().globalIndex();
            Scalar coeff;
            if(&cell == &donor)
            {
                coeff = betaFace*flux;
                centralCoeff += (1. - betaFace)*flux;
            }
            else
            {
                coeff = (1. - betaFace)*flux;
                centralCoeff += betaFace*flux;
            }

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
            eqn.boundaries()(row) -= coeff*field.prev()[nb.cell().id()];
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar flux = dot(u.faces()[bd.face().id()], bd.outwardNorm());

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= 2.*flux*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += flux;
                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("hc", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff + 2.*cell.volume()/timeStep));
        eqn.boundaries()(row) += 2.*cell.volume()*field.prev()[cell.id()]/timeStep - centralCoeff*field.prev()[cell.id()];
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, const std::vector<ImmersedBoundaryObject> &ibObjs, ScalarFiniteVolumeField &field, Scalar timeStep, Type type)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);
    VectorFiniteVolumeField gradField = grad(field);

    entries.reserve(5*field.grid.nActiveCells());

    for(const Cell &cell: field.grid.fluidCells())
    {
        Scalar centralCoeff = 0.;
        const size_t row = cell.globalIndex();

        for(const InteriorLink &nb: cell.neighbours())
        {
            //- Determine the donor and acceptor cell

            const Scalar flux = dot(u.faces()[nb.face().id()], nb.outwardNorm());
            const Cell &donor = flux > 0. ? cell : nb.cell();
            const Cell &acceptor = flux > 0. ? nb.cell() : cell;

            const Scalar gammaD = field[donor.id()];
            const Scalar gammaA = field[acceptor.id()];
            const Scalar gammaU = std::max(std::min(gammaA + 2*dot(donor.centroid() - acceptor.centroid(), gradField[donor.id()]), 1.), 0.);

            const Scalar coD = sqrt(u.faces()[nb.face().id()].magSqr()*timeStep*timeStep/nb.rCellVec().magSqr());

            const Scalar gammaTilde = (gammaD - gammaU)/(gammaA - gammaU);

            Scalar gammaTildeF;

            switch(type)
            {
            case HC:
                gammaTildeF = hc(gammaTilde, coD);
                break;
            case UQ:
                gammaTildeF = uq(gammaTilde, coD);
                break;
            }

            Scalar betaFace = (gammaTildeF - gammaTilde)/(1. - gammaTilde);

            if(isnan(betaFace)) // Upwind if a valid value is not found
                betaFace = &cell == &donor ? 0. : 1.;

            betaFace = std::max(std::min(betaFace, 1.), 0.);

            const size_t col = nb.cell().globalIndex();
            Scalar coeff;
            if(&cell == &donor)
            {
                coeff = betaFace*flux;
                centralCoeff += (1. - betaFace)*flux;
            }
            else
            {
                coeff = (1. - betaFace)*flux;
                centralCoeff += betaFace*flux;
            }

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
            eqn.boundaries()(row) -= coeff*field.prev()[nb.cell().id()];
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar flux = dot(u.faces()[bd.face().id()], bd.outwardNorm());

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= 2.*flux*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += flux;
                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("hc", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff + 2.*cell.volume()/timeStep));
        eqn.boundaries()(row) += 2.*cell.volume()*field.prev()[cell.id()]/timeStep - centralCoeff*field.prev()[cell.id()];
    }

    for(const ImmersedBoundaryObject &ibObj: ibObjs)
    {
        for(const Cell &cell: ibObj.cells())
        {
            const size_t row = cell.globalIndex();

            const Point2D &bp = ibObj.boundaryPoint(cell);
            const Point2D &ip = ibObj.imagePoint(cell);
            const Scalar gammaIp = ibObj.imagePointVal(cell, field);
            const Vector2D gradGammaIp = ibObj.imagePointVal(cell, gradField);

            const Scalar gammaBp = std::min(std::max(gammaIp + dot(gradGammaIp, bp - ip), 0.), 1.);

            const auto &kNN = ibObj.imagePointCells(cell);

            std::vector<Index> cols = {
                kNN[0].get().globalIndex(),
                kNN[1].get().globalIndex(),
                kNN[2].get().globalIndex(),
                kNN[3].get().globalIndex(),
            };

            auto coeffs = ibObj.imagePointInterpolation(cell)(ip);

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, 0.5));

            for(int i = 0; i < cols.size(); ++i)
                entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, cols[i], 0.5*coeffs[i]));

            eqn.boundaries()(row) = gammaBp;
        }
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, const VectorFiniteVolumeField &m, ScalarFiniteVolumeField &field, Scalar timeStep)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);
    VectorFiniteVolumeField gradField = grad(field);
    const Scalar k = 1.;

    entries.reserve(5*field.grid.nActiveCells());

    for(const Cell &cell: field.grid.fluidCells())
    {
        Scalar centralCoeff = 0.;
        const size_t row = cell.globalIndex();

        for(const InteriorLink &nb: cell.neighbours())
        {
            //- Determine the donor and acceptor cell

            const Scalar flux = dot(u.faces()[nb.face().id()], nb.outwardNorm());
            const Cell &donor = flux > 0. ? cell : nb.cell();
            const Cell &acceptor = flux > 0. ? nb.cell() : cell;

            const Scalar gammaD = field[donor.id()];
            const Scalar gammaA = field[acceptor.id()];
            const Scalar gammaU = std::max(std::min(gammaA + 2*dot(donor.centroid() - acceptor.centroid(), gradField[donor.id()]), 1.), 0.);

            const Scalar coD = sqrt(u.faces()[nb.face().id()].magSqr()*timeStep*timeStep/nb.rCellVec().magSqr());

            const Scalar gammaTilde = (gammaD - gammaU)/(gammaA - gammaU);

            const Scalar thetaF = acos(fabs(dot(m[cell.id()], nb.rCellVec().unitVec())));
            const Scalar psiF = std::min(k*(cos(2*thetaF) + 1.)/2., 1.);

            const Scalar gammaTildeF = psiF*hc(gammaTilde, coD) + (1. - psiF)*uq(gammaTilde, coD);

            Scalar betaFace = (gammaTildeF - gammaTilde)/(1. - gammaTilde);

            if(isnan(betaFace)) // Upwind if a valid value is not found
                betaFace = &cell == &donor ? 0. : 1.;

            betaFace = std::max(std::min(betaFace, 1.), 0.);

            const size_t col = nb.cell().globalIndex();
            Scalar coeff;
            if(&cell == &donor)
            {
                coeff = betaFace*flux;
                centralCoeff += (1. - betaFace)*flux;
            }
            else
            {
                coeff = (1. - betaFace)*flux;
                centralCoeff += betaFace*flux;
            }

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
            eqn.boundaries()(row) -= coeff*field.prev()[nb.cell().id()];
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar flux = dot(u.faces()[bd.face().id()], bd.outwardNorm());

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= 2.*flux*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += flux;
                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("hc", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff + 2.*cell.volume()/timeStep));
        eqn.boundaries()(row) += 2.*cell.volume()*field.prev()[cell.id()]/timeStep - centralCoeff*field.prev()[cell.id()];
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

} // end namespace cicsam
