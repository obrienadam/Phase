#include "ImmersedBoundary.h"
#include "Solver.h"
#include "GhostCellImmersedBoundary.h"

ImmersedBoundary::ImmersedBoundary(const Input &input, Solver &solver)
    :
      solver_(solver),
      cellStatus_(solver.addScalarField("cell_status"))
{
    try //- Lazy way to check if any immersed boundary input is present
    {
        input.boundaryInput().get_child("ImmersedBoundaries");
    }
    catch(...)
    {
        return;
    }

    for(const auto& ibObject: input.boundaryInput().get_child("ImmersedBoundaries"))
    {
        printf("Initializing immersed boundary object \"%s\".\n", ibObject.first.c_str());

        ibObjs_.push_back(
                    ImmersedBoundaryObject(ibObject.first,
                                           solver_.grid(),
                                           Vector2D(ibObject.second.get<std::string>("geometry.center")),
                                           ibObject.second.get<Scalar>("geometry.radius"))
                    );

        for(const auto& child: ibObject.second)
        {
            if(child.first == "geometry")
                continue;

            std::string type = child.second.get<std::string>("type");
            ImmersedBoundaryObject::BoundaryType boundaryType;
            Scalar boundaryRefValue = 0.;

            if(type == "fixed")
                boundaryType = ImmersedBoundaryObject::FIXED;
            else if(type == "normal_gradient")
                boundaryType = ImmersedBoundaryObject::NORMAL_GRADIENT;
            else if(type == "contact_angle")
                boundaryType = ImmersedBoundaryObject::CONTACT_ANGLE;
            else if(type == "partial_slip")
            {
                boundaryType = ImmersedBoundaryObject::PARTIAL_SLIP;
                boundaryRefValue = child.second.get<Scalar>("lambda");
            }
            else
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "unrecognized boundary type \"" + type + "\".");

            printf("Setting boundary type \"%s\" for field \"%s\".\n", type.c_str(), child.first.c_str());

            ibObjs_.back().addBoundaryType(child.first, boundaryType);
            ibObjs_.back().addBoundaryRefValue(child.first, boundaryRefValue);

            if(child.first == "p")
                ibObjs_.back().addBoundaryType("pCorr", boundaryType);
        }

        ibObjs_.back().setInternalCells();
    }

    setCellStatus();
}

Equation<ScalarFiniteVolumeField> ImmersedBoundary::eqns(ScalarFiniteVolumeField &field)
{
    return gc::ib(ibObjs_, field);
}

Equation<VectorFiniteVolumeField> ImmersedBoundary::eqns(VectorFiniteVolumeField &field)
{
    return gc::ib(ibObjs_, field);
}

//- Protected

void ImmersedBoundary::setCellStatus()
{
    for(const Cell &cell: solver_.grid().inactiveCells())
        cellStatus_[cell.id()] = SOLID;

    for(const Cell &cell: solver_.grid().fluidCells())
        cellStatus_[cell.id()] = FLUID;

    for(const ImmersedBoundaryObject& ibObj: ibObjs_)
        for(const Cell &cell: ibObj.cells())
            cellStatus_[cell.id()] = IB;
}
