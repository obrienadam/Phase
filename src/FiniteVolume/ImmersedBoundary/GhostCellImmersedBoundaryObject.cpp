#include "GhostCellImmersedBoundaryObject.h"
#include "ForcingCellStencil.h"

GhostCellImmersedBoundaryObject::GhostCellImmersedBoundaryObject(const std::string &name, Label id, FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryObject(name, id, grid)
{

}

void GhostCellImmersedBoundaryObject::update(Scalar timeStep)
{
    if(motion_)
    {
        motion_->update(*this, timeStep);
        updateCells();
    }
}

void GhostCellImmersedBoundaryObject::updateCells()
{
    fluid_->add(freshCells_);
    freshCells_.clear();
    for(const Cell& cell: cells_) //- New fresh cells
        if(!isInIb(cell.centroid()))
            freshCells_.add(cell);

    ibCells_.clear();
    solidCells_.clear();

    switch (shapePtr_->type())
    {
        case Shape2D::CIRCLE:
            for (const Cell &cell: fluid_->itemsWithin(
                    *static_cast<Circle*>(shapePtr_.get()))) //- The circle method is much more efficient
                cells_.add(cell);
            break;
        case Shape2D::BOX:
            for (const Cell &cell: fluid_->itemsWithin(
                    *static_cast<Box*>(shapePtr_.get()))) //- The box method is much more efficient
                cells_.add(cell);
            break;
        default:
            for (const Cell &cell: fluid_->itemsWithin(*shapePtr_))
                cells_.add(cell);
    }

    auto isIbCell = [this](const Cell &cell) {
        if(!isInIb(cell.centroid()))
            return false;

        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell().centroid()))
                return true;

        for (const DiagonalCellLink &dg: cell.diagonals())
            if (!isInIb(dg.cell().centroid()))
                return true;

        return false;
    };

    for (const Cell &cell: cells_)
        if (isIbCell(cell))
            ibCells_.add(cell);
        else if(isInIb(cell.centroid()))
            solidCells_.add(cell);

    grid_.setCellsActive(*fluid_);
    grid_.setCellsActive(ibCells_);
    grid_.setCellsActive(freshCells_);
    grid_.setCellsInactive(solidCells_);
    constructStencils();
}

Equation<Scalar> GhostCellImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    Equation<Scalar> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Scalar bRefValue = getBoundaryRefValue<Scalar>(field.name());

    for (const GhostCellStencil &st: stencils_)
    {
        Scalar centralCoeff;
        std::vector<Scalar> coeffs = st.ipCoeffs();

        //- Boundary assembly
        switch (bType)
        {
            case FIXED:
                centralCoeff = 0.5;
                for (Scalar &coeff: coeffs)
                    coeff *= 0.5;

                eqn.addSource(st.cell(), -bRefValue);
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                centralCoeff = -1.;
                break;

            default:
                throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }

        eqn.add(st.cell(), st.cell(), centralCoeff);

        int i = 0;
        for (const Cell &ipCell: st.ipCells())
            eqn.add(st.cell(), ipCell, coeffs[i++]);
    }

    for(const Cell& cell: freshCells_)
    {
        ForcingCellStencil st(cell, shape(), *fluid_);
        auto cells = st.nbCells();
        std::vector<Scalar> coeffs;
        Scalar src;

        switch(bType)
        {
            case FIXED:
                coeffs = st.dirichletCellCoeffs();
                src = st.dirichletBoundaryCoeff()*bRefValue;
                break;
            case NORMAL_GRADIENT:
                coeffs = st.neumannCellCoeffs();
                src = st.neumannBoundaryCoeff()*0.;
        }

        eqn.add(cell, cell, 1.);
        for(int i = 0; i < 2; ++i)
            eqn.add(cell, cells[i], -coeffs[i]);
        eqn.addSource(cell, -src);
    }

    return eqn;
}

Equation<Vector2D> GhostCellImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    Equation<Vector2D> eqn(field);
    BoundaryType bType = boundaryType(field.name());

    for (const GhostCellStencil &st: stencils_)
    {
        Scalar centralCoeff;
        std::vector<Scalar> coeffs = st.ipCoeffs();

        //- Boundary assembly
        switch (bType)
        {
            case FIXED:
                centralCoeff = 0.5;
                for (Scalar &coeff: coeffs)
                    coeff *= 0.5;

                eqn.addSource(st.cell(), -velocity());
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                centralCoeff = -1.;
                break;

            default:
                throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }

        eqn.add(st.cell(), st.cell(), centralCoeff);

        int i = 0;
        for (const Cell &ipCell: st.ipCells())
            eqn.add(st.cell(), ipCell, coeffs[i++]);
    }

    for(const Cell& cell: freshCells_)
    {
        ForcingCellStencil st(cell, shape(), *fluid_);
        auto cells = st.nbCells();
        std::vector<Scalar> coeffs;
        Vector2D src;

        switch(bType)
        {
            case FIXED:
                coeffs = st.dirichletCellCoeffs();
                src = st.dirichletBoundaryCoeff()*velocity(st.xc());
                break;
            case NORMAL_GRADIENT:
                coeffs = st.neumannCellCoeffs();
                src = st.neumannBoundaryCoeff()*Vector2D();
        }

        eqn.add(cell, cell, 1.);
        for(int i = 0; i < 2; ++i)
            eqn.add(cell, cells[i], -coeffs[i]);
        eqn.addSource(cell, -src);
    }

    return eqn;
}

void GhostCellImmersedBoundaryObject::computeNormalForce(const ScalarFiniteVolumeField &rho,
                                                         const VectorFiniteVolumeField &u,
                                                         const ScalarFiniteVolumeField &p)
{
    normalForce_ = Vector2D(0, 0);
    std::vector<std::pair<Point2D, Scalar>> forcePts;

    for(const GhostCellStencil &st: stencils())
    {
        Scalar rhoB = (st.ipValue(rho) + rho(st.cell()))/2.;
        Scalar uB = dot((st.ipValue(u) + u(st.cell()))/2., st.unitNormal());
        Scalar pB = (st.ipValue(p) + p(st.cell()))/2.;

        forcePts.push_back(std::make_pair(st.boundaryPoint(), rhoB*uB*uB/2 + pB));
    }

    //- Sort the force points
    std::sort(forcePts.begin(), forcePts.end(), [this](const std::pair<Point2D, Scalar>& lhs, const std::pair<Point2D, Scalar>& rhs){
        Vector2D rVecL = lhs.first - shapePtr_->centroid();
        Vector2D rVecR = rhs.first - shapePtr_->centroid();

        Scalar thetaL = atan2(rVecL.y, rVecL.x);
        Scalar thetaR = atan2(rVecR.y, rVecR.x);

        return (thetaL < 0 ? thetaL + 2*M_PI: thetaL) < (thetaR < 0 ? thetaR + 2*M_PI: thetaR);
    });

    //- Integrate
    for(int i = 0, end = forcePts.size(); i < end; ++i)
    {
        Point2D ptA = forcePts[i].first;
        Point2D ptB = forcePts[(i + 1)%end].first;
        Scalar fA = forcePts[i].second;
        Scalar fB = forcePts[(i + 1)%end].second;
        Point2D sb = (ptB - ptA).normalVec();

        normalForce_ += (fA + fB)/2.*sb;
    }
}

void GhostCellImmersedBoundaryObject::computeShearForce(const ScalarFiniteVolumeField &mu,
                                                        const VectorFiniteVolumeField &u)
{
    shearForce_ = Vector2D(0, 0);
}

//- Protected methods

void GhostCellImmersedBoundaryObject::constructStencils()
{
    stencils_.clear();
    for (const Cell &cell: ibCells_)
        stencils_.push_back(GhostCellStencil(cell, *shapePtr_, grid_));
}
