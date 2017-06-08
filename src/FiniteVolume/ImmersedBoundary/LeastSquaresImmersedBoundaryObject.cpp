#include <fstream>

#include "LeastSquaresImmersedBoundaryObject.h"

LeastSquaresImmersedBoundaryObject::LeastSquaresImmersedBoundaryObject(const std::string &name, Label id, FiniteVolumeGrid2D &grid)
    :
      ImmersedBoundaryObject(name, id, grid)
{

}

void LeastSquaresImmersedBoundaryObject::update(Scalar timeStep)
{
    //- Doesn't move, stationary
}

void LeastSquaresImmersedBoundaryObject::updateCells()
{
    if(updateCellsCalled_)
        throw Exception("LeastSquaresImmersedBoundaryObject", "updateCells", "cannot call updateCells more than once.");

    ibCells_.clear();
    solidCells_.clear();

    //- Find the solid cells
    switch (shapePtr_->type())
    {
    case Shape2D::CIRCLE:
        for (const Cell &cell: fluid_->itemsWithin(
                 *(Circle *) shapePtr_.get())) //- The circle method is much more efficient
            cells_.add(cell);
        break;
    case Shape2D::BOX:
        for (const Cell &cell: fluid_->itemsWithin(
                 *(Box *) shapePtr_.get())) //- The box method is much more efficient
            cells_.add(cell);
        break;
    default:
        for (const Cell &cell: fluid_->itemsWithin(*shapePtr_))
            cells_.add(cell);
    }

    solidCells_.add(cells_);

    //- Find the ib cells
    for(const Cell& cell: solidCells_)
        for(const InteriorLink& nb: cell.neighbours())
            if(!shape().isInside(nb.cell().centroid()))
            {
                cells_.add(nb.cell());
                ibCells_.add(nb.cell());
            }

    grid_.setCellsActive(ibCells_);
    grid_.setCellsInactive(solidCells_);
    constructStencils();
    updateCellsCalled_ = true;
}

Equation<Scalar> LeastSquaresImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    Equation<Scalar> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Scalar bRefValue = getBoundaryRefValue<Scalar>(field.name());

    for (const LeastSquaresStencil &st: stencils_)
    {
        //- Boundary assembly
        switch (bType)
        {
        case FIXED:
        {
            std::vector<Scalar> coeffs = st.dirichletCellCoeffs();
            std::vector<Ref<const Cell>> cells = st.iCells();

            eqn.add(st.cell(), st.cell(), 1.);

            for(int i = 0; i < cells.size(); ++i)
                eqn.add(st.cell(), cells[i], -coeffs[i]);

            coeffs = st.dirichletBoundaryCoeffs();

            for(int i = 0; i < coeffs.size(); ++i)
                eqn.addSource(st.cell(), -coeffs[i]*bRefValue);
        }

            break;

        case NORMAL_GRADIENT:
        {
            std::vector<Scalar> coeffs = st.neumannCellCoeffs();
            std::vector<Ref<const Cell>> cells = st.iCells();

            eqn.add(st.cell(), st.cell(), 1.);

            for(int i = 0; i < cells.size(); ++i)
                eqn.add(st.cell(), cells[i], -coeffs[i]);

            coeffs = st.neumannBoundaryCoeffs();

            for(int i = 0; i < coeffs.size(); ++i)
                eqn.addSource(st.cell(), -coeffs[i]*bRefValue);
        }
            break;

        default:
            throw Exception("LeastSquaresImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }
    }

    return eqn;
}

Equation<Vector2D> LeastSquaresImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    Equation<Vector2D> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Vector2D bRefValue = getBoundaryRefValue<Vector2D>(field.name());

    for (const LeastSquaresStencil &st: stencils_)
    {
        //- Boundary assembly
        switch (bType)
        {
        case FIXED:
        {
            std::vector<Scalar> coeffs = st.dirichletCellCoeffs();
            std::vector<Ref<const Cell>> cells = st.iCells();

            eqn.add(st.cell(), st.cell(), 1.);

            for(int i = 0; i < cells.size(); ++i)
                eqn.add(st.cell(), cells[i], -coeffs[i]);

            coeffs = st.dirichletBoundaryCoeffs();

            for(int i = 0; i < coeffs.size(); ++i)
                eqn.addSource(st.cell(), -coeffs[i]*bRefValue);
        }

            break;

        case NORMAL_GRADIENT:
        {
            std::vector<Scalar> coeffs = st.neumannCellCoeffs();
            std::vector<Ref<const Cell>> cells = st.iCells();

            eqn.add(st.cell(), st.cell(), 1.);

            for(int i = 0; i < cells.size(); ++i)
                eqn.add(st.cell(), cells[i], -coeffs[i]);

            coeffs = st.neumannBoundaryCoeffs();

            for(int i = 0; i < coeffs.size(); ++i)
                eqn.addSource(st.cell(), -coeffs[i]*bRefValue);
        }
            break;

        default:
            throw Exception("LeastSquaresImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }
    }

    return eqn;
}

void LeastSquaresImmersedBoundaryObject::constructStencils()
{
    for(const Cell& cell: ibCells_)
        stencils_.push_back(LeastSquaresStencil(cell, *shapePtr_));

    std::ofstream fout("stencils.dat");

    for(const LeastSquaresStencil& st: stencils_)
    {
        fout << "Geometry x=0 y=0 T=LINE C=BLACK LT=0.2 CS=GRID\n"
             << st.iCells().size() + st.boundaryPoints().size() << "\n";

        for(const Cell& cell: st.iCells())
        {
            fout << "2\n"
                 << st.cell().centroid().x << " " << st.cell().centroid().y << "\n"
                 << cell.centroid().x << " " << cell.centroid().y << "\n";
        }

        for(const auto& bp: st.boundaryPoints())
        {
            fout << "2\n"
                 << st.cell().centroid().x << " " << st.cell().centroid().y << "\n"
                 << bp.first.x << " " << bp.first.y << "\n";
        }
    }

    fout.close();
}
